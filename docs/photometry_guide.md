# Pipeline — Technical Notes

---

## Version history

| Version | Changes |
|---------|---------|
| v0.1.7 | Photometry on the Bayer green channel (`seqextract_Green` before registration, AAVSO filter `TG`); plate solve of the reference frame only; green-channel gain; FWHM in arcsec from the plate solution |
| v0.1.6 | AAVSO export accepted by WebObs (commented column line, `STD` apparent magnitudes, observer code); V_app from the flux mean of the comp stars Siril used |
| v0.1.5 | PRIMARY path via `findcompstars` + `-ninastars`; ensemble V_app; FWHM from `.seq` |
| v0.1.4 | Frame border safety margin raised to 200 px |
| v0.1.3 | In-frame filtering of comp stars for FALLBACK A/B |
| v0.1.2 | Early frame check at step 5; removed broken Siril annotate |
| v0.1.1 | WCS-based frame filtering |
| v0.1.0 | Initial release |

---

## v0.1.7 — How it works

### 1. Session path

The user provides the session folder path, which must contain a `lights/` subfolder with Seestar FITS frames. On load, the script reads the FITS headers of the first file to get:

- `OBJCTRA` / `OBJCTDEC` or `RA` / `DEC` → field center (equatorial coordinates)
- `NAXIS1` / `NAXIS2` + plate scale → field size in degrees
- Count of `.fit` files in `lights/`

### 2. VSX catalog query

A cone search is sent to **VizieR** (`B/vsx/vsx`) around the field center. Search radius is half the field diagonal, rounded up to the nearest 0.1° — wide enough to catch stars near the edges.

Default magnitude limit: **14**. All VSX variable stars within the cone and below this threshold are returned and displayed in the table. The user can filter by name or type in the UI and adjust the magnitude limit.

### 3. Target and comparison star selection

After selecting a target from the table, the script uses the comparison and check stars chosen in Argos when the session has an Argos selection, and otherwise calls `findcompstars` (Siril) or queries VizieR APASS DR9 depending on availability. See [photometry path selection](#photometry-path-selection) for details.

### 4. Filter and calibration frames

**Filter:** the photometry is done on the Bayer green channel, so the AAVSO `FILT` field is `TG` whatever the filter wheel position. Choosing the LP option only adds a note to the export.

**Obscode:** your AAVSO observer code. Without it no `aavso.csv` is written (see [AAVSO export](#aavso-export)).

**Mount — Equatorial mode:** tick it for a session shot in equatorial mode. In alt-az mode the field rotates during the session, so only the inscribed circle of the frame (radius `min(W, H)/2 − 50 px`) is covered by every registered frame: the VSX list keeps only stars inside it, and step 5 stops if the target is outside. In equatorial mode there is no field rotation, so the whole frame is usable: the VSX list keeps stars at least 50 px from the edge, and step 5 only needs the target's sky annulus inside the frame (`max(35, outer + 5)` px from the edge). The setting is not saved: tick it again when you reopen the session, then re-query VSX.

**Calibration:** the user can provide folders of darks, flats, and/or bias frames. If at least one type is given, the pipeline creates masters automatically at step 0:

```
link dark_ → stack dark_ rej 3 3 -nonorm  → master_dark.fit
link flat_  → stack flat_  rej 3 3 -norm=mul → master_flat.fit
link bias_  → stack bias_  rej 3 3 -nonorm  → master_bias.fit
calibrate light_ -dark=... -flat=... -bias=... -cc=dark -cfa  → pp_light_*.fit
```

Stacking uses Winsorized sigma clipping 3σ, robust against cosmic rays and satellite trails. Without calibration, the active sequence stays `light_`.

### 5. Run — steps and algorithms

**Step 1 — Sequence and green channel**

```
link light -out=process/
seqextract_Green light_        →  Green_light_*.fit   (or Green_pp_light_ after calibration)
```

`seqextract_Green` keeps the two green pixels of each GRBG cell and averages them, (G1+G2)/2, in a half-resolution image (1080 × 1920 px, 5.8 µm pixels, 7.35″/px). The photometric band is the Bayer green (AAVSO `TG`), the same plane as the Argos live curve. `BAYERPAT` is removed and `DATE-OBS`, `EXPTIME`, `AIRMASS` are kept.

Before v0.1.7 the frames stayed CFA: registration interpolated the mosaic as a mono image and the aperture summed red, green and blue pixels — a broad band reported as `CV`.

**Steps 2+3 — Registration**

```
register Green_light_ -2pass
seqapplyreg Green_light_ -framing=min -filter-round=2.5k  →  r_Green_light_*.fit
```

`register -2pass` detects the stars, picks the reference frame (FWHM, star count) and computes a homography per frame, written to the `.seq` only. `seqapplyreg` writes the aligned frames (Lanczos4 interpolation), cropped to the area common to all frames (`-framing=min`); `-filter-round=2.5k` drops frames whose roundness is more than 2.5σ below the median.

**Step 4 — Plate solve of the reference frame**

```
load r_Green_light_<ref>
platesolve -force -noflip -focal=160 -pixelsize=5.8
save r_Green_light_<ref>
```

The frames are aligned on the reference frame, so its WCS places the stars on all of them: `light_curve` projects the `-ninastars` coordinates on the reference image. `seqplatesolve` is not used: it drops from the sequence every frame it cannot solve — about half of the undersampled green frames in a test session.

**Step 5 — Photometry**

```
setphot -aperture=10 -inner=10 -outer=15 -dyn_ratio=4 -min_val=-1000 -max_val=60000 -gain=<2 × EGAIN>
light_curve r_Green_light_ 0 -ninastars=comp_stars.csv
```

The aperture is dynamic, `0.5 × dyn_ratio × FWHM` = 2 FWHM for each star on each frame (`-aperture` only applies to a forced radius). The sky annulus is 10–15 green pixels (73″–110″). `light_curve` fits each star, measures it in the aperture and computes V−C against the mean flux of the comparison stars. Raw output is `light_curve.dat` (JD at mid-exposure, V−C, error).

**Python post-processing**

After Siril exits, the Python pipeline:
1. Reads `light_curve.dat`
2. Computes `V_app = V_C + C_cat`, `C_cat` being the flux mean of the catalogue V magnitudes of the comp stars Siril used (see [V_app](#v_app-and-the-lp-filter-bias))
3. Extracts per-frame FWHM from the `.seq` R0 lines (mean and weighted FWHM, in pixels), aligns with `DATE-OBS` from FITS headers, converts to arcsec with the scale of the plate solution
4. Writes result files

### 6. Output

```
results/StarName/
├── photometry.csv       JD, V_C, V_app, err
├── fwhm.csv             JD, FWHM_arcsec, wFWHM_arcsec
└── aavso.csv            AAVSO Extended Format (only when it can be submitted)
```

`photometry.csv`:
```
# Ensemble V (flux mean of the APASS comp stars): 12.284
# V_app = V_C + ensemble_V
JD,V_C,V_app,err
2461161.334643,1.2421,13.5261,0.0868
```

`fwhm.csv`:
```
# Plate scale: 7.3488 arcsec/px
JD,FWHM_arcsec,wFWHM_arcsec
2461309.336766,11.896,15.215
```

---

## Photometry path selection

**ARGOS** — the stars chosen in Argos + `-ninastars`

Argos writes the comparison and check stars it used during acquisition, taken from the target's AAVSO VSP sequence, to `photometry_selection.json` in its session folder. The script fills the **Argos:** field with that file when it loads the session; browse to another one, or clear the field to use the automatic paths below.

The comparison stars go to Siril as `Comp2` (AAVSO) rows of `comp_stars.csv`, named by AUID. Only those inside the reference frame and with a catalogue V magnitude are kept. The check star is measured in a second `light_curve` run against the same comparison stars, since `light_curve` only reports its target; `aavso.csv` gets its magnitude in `KMAG` only if Siril used the same comparison stars in both runs.

The selection is not used — and the log says why — when its target is more than 30″ from the selected star, when fewer than two comparison stars remain, or when the reference frame has no WCS. The script then falls back to PRIMARY.

**PRIMARY** — `findcompstars` + `-ninastars`

Siril handles the sky-to-pixel coordinate conversion internally via the WCS headers. Requires valid WCS on all frames and at least 3 comp stars in the field.

**FALLBACK A** — pixel coordinates from VizieR

Comp stars queried from VizieR APASS DR9 (`e_Vmag < 0.05`, `|delta_V| < 2.0`), converted to pixels via TAN projection from the reference frame WCS. Coordinates rounded to integers — Siril 1.4.x rejects decimals.

**FALLBACK B** — single nearest bright star

Last resort. Differential magnitude is noisier, and any intrinsic variability of the comp star goes undetected.

---

## V_app and the LP filter bias

`V_app = V_C + C_cat`, with `C_cat = -2.5 log10(mean(10^(-0.4 V_i)))` over the catalogue V magnitudes of the comp stars.

Siril's `light_curve` computes `V_C` against the mean *flux* of the comp stars (`photometry.c`, `new_light_curve`), so the catalogue magnitudes must be averaged the same way. A median of the V magnitudes is biased when the comp stars span a few magnitudes — about 0.6 mag for comps between V=10.2 and 13.1, the brightest dominating the flux mean (earlier versions had this bias).

`C_cat` must cover exactly the stars Siril averaged. Siril skips `-ninastars` stars on the frame border (logged by name) and drops any comp star measured on fewer than 4/5 of the frames (only the count is logged). The pipeline reads both from Siril's output; when a star was dropped it cannot tell which, so V_app — and `aavso.csv` — are left out, and the log says so. Siril's output is forced to English for this (`LANGUAGE=C`); a language chosen in Siril's own preferences overrides that, and only English and French messages are recognised.

V_C is unaffected by the LP filter (target and comp stars go through the same filter). The bias enters via `C_cat`, which is in standard APASS V while the observation is in LP. For a cataclysmic variable (blue + red) vs. G/K comp stars, expect ±0.2–0.5 mag absolute offset. Relative variations within a session are reliable.

---

## AAVSO export

`aavso.csv` follows the [AAVSO Extended File Format](https://www.aavso.org/aavso-extended-file-format), to upload on WebObs. It is written only when it can be submitted as is: an observer code is set, V_app is available (AAVSO cannot use the differential V_C) and the star name fits the format. Otherwise the log says why and any previous `aavso.csv` is removed.

```
#TYPE=EXTENDED
#OBSCODE=ABC
#SOFTWARE=Siril+seestar_varstar_siril.py v0.1.7
#DELIM=,
#DATE=JD
#OBSTYPE=CCD
#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
ES UMa,2461161.334643,11.707,0.021,TG,NO,STD,ENSEMBLE,na,na,na,1.302,na,APASS DR9,Seestar S30 Pro; Siril ensemble of 3 APASS DR9 comparison stars (flux-mean V=12.119); no check star; LP_filter_Seestar_S30Pro
```

| Field | Value | Why |
|-------|-------|-----|
| column-name line | starts with `#` | WebObs reads every other line as an observation |
| NAME | VSX name without `V* ` | AAVSO star identifier, max 30 characters |
| DATE | JD UTC, mid-exposure | Siril adds half the exposure time |
| MAG / MERR | V_app / error, 3 decimals | NaN and MERR > 0.5 excluded |
| TRANS / MTYPE | `NO` / `STD` | untransformed, standardised on the comp stars' catalogue magnitudes |
| CNAME / CMAG | `ENSEMBLE` / `na` | the spec's ensemble convention |
| KNAME / KMAG | Argos check star AUID / its magnitude, else `na` | the spec recommends a check star for ensemble photometry; the automatic paths have none |
| GROUP | `na` | time series, single filter |
| CHART | VSP chart ID (Argos stars), else `APASS DR9` | a VSP sequence is reported by its chart ID; self-chosen comp stars by their catalogue |
| NOTES | comp source and count, `C_cat`, check star, filter note | the spec asks for information on self-chosen comp stars |

`tests/aavso_spec.py` checks exported files against the rules of the spec, and is itself checked against the spec's example files.

`seqsetmag` (Siril's magnitude calibration) is not scriptable in headless mode — V_app is therefore computed externally.

---

## Siril 1.4.x bugs worked around

| Bug / limitation | Workaround |
|-----------------|------------|
| `-at` pixel coords must be integers | `round()` before building the command |
| `-autoring` incompatible with `-at` | Fixed `setphot` values instead |
| `seqsetmag` not scriptable in headless | V_app computed from `comp_stars.csv` |
| `.ssf` parser includes quotes in path | Always use `-out=path` unquoted |
| `seqpsf` headless: console output only, no file | Not used; FWHM read from `.seq` instead |
| `light_curve -wcs` can fail on frames with missing WCS | PRIMARY uses `-ninastars`; FALLBACK A uses `-at` |

---

## Plate scale

```
206.265 * 2.9 µm / 160 mm = 3.74 arcsec/px   (raw Bayer pixel, nominal)
206.265 * 5.8 µm / 160 mm = 7.48 arcsec/px   (green-channel pixel, nominal)
```

The plate solution of a test session gives 3.673″/px on the raw mosaic and 7.349″/px on the green channel (effective focal length ≈ 163 mm). `fwhm.csv` uses the scale of the reference frame's WCS.

## Gain e⁻/ADU

`get_gain_eadu()` tries headers in order: `EGAIN`, `EPERDN`, `GAIN_E`, `CCDGAIN`, then `GAIN`. Value accepted only if `0.05 < g < 30`; otherwise 1.0 e⁻/ADU. Argos writes `EGAIN` (the Alpaca driver value, 0.8 e⁻/ADU at gain 80).

A green-channel pixel is the mean of two raw pixels: one of its ADU holds twice the electrons, so `setphot -gain` receives twice the raw gain. `seqextract_Green` drops `EGAIN`, so Siril uses this value.
