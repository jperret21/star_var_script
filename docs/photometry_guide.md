# Variable Star Photometry — Technical Reference

Science and implementation details behind `seestar_varstar_siril.py`.

---

## 1. Differential photometry

Differential photometry measures the brightness of a target star **relative to comparison stars** in the same image. Because all stars are affected identically by atmospheric extinction, thin clouds, and varying seeing within a single frame, the differential magnitude is stable even under imperfect conditions:

```
diff_mag(target) = −2.5 × log10(flux_target / flux_ensemble)
```

where `flux_ensemble` is a weighted combination of all comparison stars.

**What it cancels:**
- First-order atmospheric extinction (same air mass for all stars)
- Variable transparency (clouds, haze)
- Detector gain variations frame-to-frame

**What it does NOT cancel:**
- Color-dependent extinction (differential color extinction) — minimized by choosing comparison stars with similar B−V color to the target
- Instrumental color response differences between the target and comparison stars (this is what a photometric transformation corrects)

The pipeline uses Siril's ensemble method: `CNAME=ENSEMBLE`, `MTYPE=DIFF`.

---

## 2. The LP filter and photometric bands

The Seestar S30 Pro ships with a built-in **light pollution (LP) filter** that blocks specific artificial lighting wavelengths (Na 589 nm, Hg 436/546 nm, etc.) while passing most of the visible spectrum.

### Why this matters for photometry

Standard photometric catalogs (APASS, AAVSO) publish magnitudes in defined bands:

| Band | λ_eff | Width | Notes |
|------|-------|-------|-------|
| Johnson B | 440 nm | 90 nm | Blue |
| Johnson V | 550 nm | 90 nm | Green-yellow, most common |
| Johnson R | 640 nm | 150 nm | Red |
| Sloan r' | 620 nm | 140 nm | Similar to R |

The LP filter's passband does not match any of these. Comparison star magnitudes from APASS are in V (or B, r, i) — applying them directly to LP-filtered measurements introduces a systematic offset that varies with the color of each star.

### AAVSO filter code: `CV`

The closest AAVSO-recognized code for LP-filtered or unfiltered OSC observations is `CV` (Clear/Visual). This is what the pipeline uses, with a note in the NOTES field:

```
NOTES = seestar_s30pro|LP_filter_Seestar_S30Pro
```

**Consequence for submitted data:** AAVSO accepts `CV` observations. They are not directly comparable to V-band data without a color transformation coefficient (Tv), which requires observations of photometric standard fields. For variability monitoring (period, amplitude, event timing), this is generally sufficient.

### Removing the filter

For future runs without the LP filter, select **"Clear — no filter (CV)"** in the Filter dropdown. The AAVSO code remains `CV` but the note is omitted, which is accurate for an unfiltered OSC observation.

---

## 3. Gain and photometric error bars

Siril's photometric error model combines:

```
σ_total² = σ_photon² + σ_sky² + σ_readout²

σ_photon = √(flux / gain)           [photon shot noise]
σ_sky     = A × σ_bg / √A_annulus   [background estimation]
σ_readout = readnoise / gain         [detector electronics]
```

where `gain` is in electrons per ADU.

### The GAIN=200 problem

Seestar FITS headers contain `GAIN=200`. This is **not** the gain in e⁻/ADU — it is the camera gain setting (analogous to ISO 200). The Sony IMX585 at the Seestar's default mode runs at approximately **0.5–1.0 e⁻/ADU**.

If Siril receives `gain=200`, it assumes 200 electrons are collected per ADU, making shot noise appear 14× smaller than it really is → **error bars are significantly underestimated**.

The pipeline reads `EGAIN` from the FITS header (the actual e⁻/ADU value, written by some processing software). If absent, it falls back to **1.0 e⁻/ADU**, which is conservative and realistic for the IMX585. This is logged at the start of every run:

```
Gain: 1.0 e-/ADU (from FITS EGAIN header, fallback=1.0)
```

For better accuracy, run the Seestar data through a tool that writes `EGAIN` to the header (Siril's calibration step does this when it knows the sensor parameters).

---

## 4. Channel selection: why channel 0

After `register -2pass` + `seqapplyreg` on raw CFA (Bayer) frames, Siril outputs **1-layer mono FITS files**. The Bayer mosaic is preserved as a single-channel image (the color information is encoded spatially, not in separate planes).

Siril's `light_curve` command uses a 0-based channel index:
- `0` → the only channel in a 1-layer sequence ✓
- `1` → would require a 3-channel RGB sequence ✗

The pipeline uses `light_curve r_light_ 0`.

### Why not debayer first?

For aperture photometry, debayering is not necessary. The key requirement is that pixel values are proportional to flux. On a CFA image, each aperture still contains many pixels of each color; the PSF-weighted average gives a good flux estimate. Debayering would:
- Add interpolation artifacts near stars
- Expand the PSF slightly
- Require more processing time

For best accuracy with a specific band (V, R), use a real photometric filter instead of debayering.

---

## 5. Pipeline step by step

### Step 0 — Calibration (optional)

If darks, flats, or bias frames are provided:

```
# For each calibration type (dark / flat / bias):
link dark -out=/path/to/masters      # sequence named dark_
stack dark_ rej 3 3 -nonorm          # Winsorized sigma clipping
  → master_dark.fit

# Apply to lights:
calibrate light_ -dark=masters/master_dark -flat=masters/master_flat
  → pp_light_*.fit
```

Stacking uses 3-σ rejection (Winsorized), which is robust against cosmic rays and satellite trails. Flats are normalized multiplicatively (`-norm=mul`).

After calibration the active sequence is `pp_light_`; without calibration it stays `light_`.

### Steps 1+2 — Link + Register (same siril-cli session)

```
link light -out=/path/to/process     # copies + names sequentially
register light_ -2pass               # star-matching, no internet needed
```

`register -2pass` works offline using detected star patterns. It computes a similarity transform (translation + rotation + scale) for each frame relative to a reference frame. Results are stored in `light_.seq`.

**Why both commands run in the same siril-cli session:** `link` creates the sequence in Siril's in-memory state. A second `siril-cli` invocation would not find the sequence file if `link` created it in an unexpected location. Running them together avoids this.

> **Known quirk:** Siril's `.ssf` parser does NOT strip quotes from parameter values like `-out="path"`. Always use `-out=path` (unquoted) for paths without spaces.

### Step 3 — seqapplyreg

```
seqapplyreg light_ -framing=max -filter-round=2.5k
  → r_light_*.fit
```

- `-framing=max`: keeps the maximum sky area (union of all frames)
- `-filter-round=2.5k`: rejects the worst 2500 frames by star elongation — discards frames with tracking errors, wind shake, or satellite trails

Output frames are interpolated to a common pixel grid. Pixel values are resampled with Lanczos interpolation (default in Siril).

### Step 4 — seqplatesolve

```
seqplatesolve r_light_ -nocache -force -focal=160 -pixelsize=2.9 -radius=2.5
  → WCS headers written to each r_light_*.fit
```

Each frame gets an independent plate solve against Gaia DR3 (online by default). The solution writes `CRVAL1/2`, `CRPIX1/2`, `CD1_1` etc. into the FITS header.

- `-nocache`: don't reuse cached solutions from a previous run
- `-force`: re-solve even if WCS already present
- `-radius=2.5`: search radius in degrees around the nominal pointing

**Why per-frame WCS?** After `seqapplyreg`, frames are aligned but the WCS of the reference frame doesn't apply to the others without individual solutions (because `-framing=max` may add padding). Per-frame WCS allows `light_curve` to find the target and comparison stars by sky coordinates in every frame.

### Step 5 — setphot + light_curve

```
setphot -aperture=8 -inner=14 -outer=21 -gain=1.0
light_curve r_light_ 0 -at=x,y -refat=x1,y1 ...
```

The pipeline uses pixel coordinates (`-at/-refat`) rather than sky coordinates (`-wcs/-refwcs`). After `seqapplyreg`, all frames share the same pixel grid, so fixed pixel positions are correct for every frame. The pixel coordinates are computed from the reference frame's WCS using a TAN gnomonic projection (SIP distortions are ignored; corrections are sub-pixel at the Seestar's 3.74 arcsec/pixel scale). This approach also bypasses Siril's internal WCS availability check, which can fail when frames have WCS from a non-Siril solver.

**Siril 1.4.3 known constraints with `light_curve -at`:**
- Pixel coordinates must be **integers** — decimal values cause `invalid arguments`.
- `-autoring` is **incompatible** with `-at` mode (Siril bug: autoring's Findstar validation rejects pixel coordinates, reporting "coordinates not in image"). Use fixed `setphot` aperture values calibrated to the expected FWHM instead (aperture ≈ 2.5 × FWHM, inner ≈ 4.2 × FWHM, outer ≈ 6.3 × FWHM).

`setphot` configures the aperture geometry (in pixels):

```
[  outer annulus 25px  ]
[  inner annulus 15px  ]
[    aperture 10px     ]   ← star flux summed here
[    star PSF          ]
```

At the Seestar's 3.74 arcsec/pixel scale:
- Aperture: 10 px = 37.4 arcsec radius
- Sky annulus: 15–25 px = 56–94 arcsec

With `-autoring`, Siril overrides these with values derived from the measured FWHM of stars in the frame:
- aperture ≈ 2.5 × FWHM
- inner ≈ 4.2 × FWHM
- outer ≈ 6.3 × FWHM

This is generally better than fixed values and adapts to varying seeing.

---

## 6. Selecting comparison stars

The pipeline fetches comparison stars from **VizieR APASS DR9** (catalog `II/336/apass9`).

Filters applied automatically:
- `e_Vmag < 0.05` — only stars with < 0.05 mag catalog uncertainty
- `|Vmag − target_mag| < 2.0` — within 2 mag of the target
- Separation from target > 5 arcsec — avoids blending

The script selects the 6 closest-in-magnitude stars. For most targets this gives a good ensemble. If fewer than 3 are found, the pipeline warns but continues.

### Manual selection

If APASS gives poor results (sparse field, crowded region, bright target):
1. Use [AAVSO VSP](https://www.aavso.org/apps/vsp) to get a chart with labeled comp stars
2. The GUI currently shows APASS auto-selection only — manual comparison star input is a planned feature

---

## 7. AAVSO Extended Format — field reference

```csv
#TYPE=EXTENDED
#OBSCODE=XXXX          ← replace with your real code
#SOFTWARE=Siril+seestar_varstar_siril.py
#FILTER=CV
#DELIM=,
#DATE=JD
#OBSTYPE=CCD
NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
SN 1993J,2460796.410700,12.3420,0.0210,CV,NO,DIFF,ENSEMBLE,na,na,na,na,1,na,seestar_s30pro|LP_filter_Seestar_S30Pro
```

| Field | Value | Rationale |
|-------|-------|-----------|
| `FILT` | `CV` | LP filter is not a standard band; CV is the correct AAVSO code for broadband/unfiltered |
| `TRANS` | `NO` | No color transformation coefficient applied |
| `MTYPE` | `DIFF` | Differential — magnitudes are relative to comparison ensemble, not transformed to standard system |
| `CNAME` | `ENSEMBLE` | Multiple comparison stars used |
| `CMAG` | `na` | Ensemble magnitude not reported as single value |
| `AMASS` | `na` | Airmass not computed (requires observer coordinates) |
| `NOTES` | `LP_filter_Seestar_S30Pro` | Documents the filter for data quality assessment |

---

## 8. Siril commands reference (1.4 syntax)

```
link basename -out=directory
    Create a numbered sequence from all FITS in current directory.
    NOTE: -out= value must NOT be quoted in .ssf scripts —
    Siril's parser includes the quotes literally in the path.

register sequence -2pass
    Star-pattern matching registration. Offline, robust.
    Writes transformation matrices to .seq file.

seqapplyreg sequence -framing=max -filter-round=Nk
    Apply registration transforms. -framing=max uses maximum field.
    -filter-round=2.5k keeps best 2500 frames by star roundness.

calibrate sequence [-bias=file] [-dark=file] [-flat=file] [-cc=banding]
    Apply calibration frames. -cc=banding corrects column banding (CMOS).
    File paths are relative to CWD, no extension.

stack sequence [type] [lo] [hi] [-nonorm | -norm=mul] [-out=name]
    Stack a sequence. rej 3 3 = Winsorized sigma clipping 3σ.
    -nonorm for darks/bias, -norm=mul for flats.

seqplatesolve sequence -nocache -force -focal=mm -pixelsize=um -radius=deg
    Per-frame plate solve using Gaia DR3 (online).
    Writes full WCS to each FITS header.

setphot -aperture=px -inner=px -outer=px -gain=value
    Configure aperture photometry. gain in e-/ADU.
    Called once before light_curve.

light_curve sequence channel [-autoring] -wcs=RA,Dec [-refwcs=RA,Dec ...]
    Aperture photometry on every frame.
    channel: 0 for 1-layer (CFA or mono), 1/2 for 3-channel RGB.
    -autoring: derive aperture radii from frame FWHM (recommended).
    -wcs: target sky coordinates (decimal degrees).
    -refwcs: comparison star coordinates (repeat per star).
    Output: light_curve.dat + light_curve.png in CWD.
```

Full reference: [siril.readthedocs.io/en/latest/Commands.html](https://siril.readthedocs.io/en/latest/Commands.html)

---

## 9. Typical precision and observing strategy

### Expected precision with Seestar S30 Pro

| Target magnitude | Expected σ (20s subs) | Notes |
|-----------------|----------------------|-------|
| < 10 mag | > 0.05 mag | Near saturation — shorten exposure |
| 10–13 mag | 0.01–0.03 mag | Good regime for Seestar |
| 13–14 mag | 0.03–0.07 mag | Usable for bright variables |
| > 14 mag | > 0.1 mag | Marginal — use calibration frames + longer subs |

### Short-period variables (RR Lyrae, δ Scuti, eclipsing binaries)

- Use raw 20s frames — don't stack
- Cover at least 1–2 full periods
- Cadence: one point every ~25s (frame + readout)

### Long-period variables (Mira, SR, semi-regular, T > 10 days)

- Stack groups of 5–10 frames to reduce scatter
- 1–3 measurements per night is sufficient
- Use **"Photometry only"** on repeat nights after the first full run

### SN monitoring (e.g. SN 1993J in M81)

- Individual frames give one point per 20s — cadence is excellent
- Target is ~12 mag in M81, well within Seestar range
- Watch for slow magnitude change over days/weeks by comparing sessions
