# Variable Star Photometry — Technical Reference

Science and implementation details behind `seestar_varstar_siril.py`.

---

## Contents

1. [Differential photometry](#1-differential-photometry)
2. [The LP filter and photometric bands](#2-the-lp-filter-and-photometric-bands)
3. [Gain and photometric error bars](#3-gain-and-photometric-error-bars)
4. [Channel selection: why channel 0](#4-channel-selection-why-channel-0)
5. [Pipeline step by step](#5-pipeline-step-by-step)
6. [Selecting comparison stars](#6-selecting-comparison-stars)
7. [Apparent magnitude computation](#7-apparent-magnitude-computation)
8. [FWHM extraction](#8-fwhm-extraction)
9. [AAVSO Extended Format — field reference](#9-aavso-extended-format--field-reference)
10. [Post-processing: varstar_postprod](#10-post-processing-varstar_postprod)
11. [Siril commands reference](#11-siril-commands-reference-14-syntax)
12. [Typical precision and observing strategy](#12-typical-precision-and-observing-strategy)
13. [Troubleshooting](#13-troubleshooting)

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
- Instrumental color response differences (this is what a photometric transformation corrects)

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

The LP filter's passband does not match any of these. Comparison star magnitudes from APASS are in V (or B, r, i) — applying them directly to LP-filtered measurements introduces a systematic offset that varies with the color of each star. For a cataclysmic variable (blue continuum + red companion) vs. G/K comp stars, this is typically ±0.2–0.5 mag.

### AAVSO filter code: `CV`

The closest AAVSO-recognized code for LP-filtered or unfiltered OSC observations is `CV` (Clear/Visual):

```
NOTES = seestar_s30pro|LP_filter_Seestar_S30Pro
```

**Consequence for submitted data:** AAVSO accepts `CV` observations. They are not directly comparable to V-band data without a color transformation coefficient (Tv). For variability monitoring (period, amplitude, event timing), this is generally sufficient.

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

The pipeline reads `EGAIN` from the FITS header (the actual e⁻/ADU value). If absent, it falls back to **1.0 e⁻/ADU**, which is conservative and realistic for the IMX585.

---

## 4. Channel selection: why channel 0

After `register -2pass` + `seqapplyreg` on raw CFA (Bayer) frames, Siril outputs **1-layer mono FITS files**. The Bayer mosaic is preserved as a single-channel image.

Siril's `light_curve` command uses a 0-based channel index:
- `0` → the only channel in a 1-layer sequence ✓
- `1` → would require a 3-channel RGB sequence ✗

The pipeline uses `light_curve r_light_ 0`.

---

## 5. Pipeline step by step

### Step 0 — Calibration (optional)

```
link dark -out=masters      → dark_
stack dark_ rej 3 3 -nonorm → master_dark.fit

calibrate light_ -dark=masters/master_dark -flat=masters/master_flat
  → pp_light_*.fit
```

Stacking uses 3-σ Winsorized sigma clipping. After calibration the active sequence is `pp_light_`; without it stays `light_`.

### Steps 1+2 — Link + Register

```
link light -out=/path/to/process
register light_ -2pass
```

`register -2pass` works offline using star-pattern matching. Writes transformation matrices to `light_.seq`.

### Step 3 — seqapplyreg

```
seqapplyreg light_ -framing=max -filter-round=2.5k  →  r_light_*.fit
```

- `-framing=max`: maximum sky area (union of all frames)
- `-filter-round=2.5k`: rejects worst frames by star elongation

### Step 4 — seqplatesolve

```
seqplatesolve r_light_ -nocache -force -focal=160 -pixelsize=2.9 -radius=2.5
```

Each frame gets an independent plate solve against Gaia DR3. Writes full WCS headers (`CRVAL`, `CRPIX`, `CD` matrix) to every `r_light_*.fit`.

**Why per-frame WCS?** After `seqapplyreg` with `-framing=max`, the reference frame WCS doesn't apply to the others without individual solutions.

### Step 5 — setphot + light_curve

```
setphot -aperture=8 -inner=14 -outer=21 -gain=1.0
light_curve r_light_ 0 -ninastars=comp_stars.csv
```

**Siril 1.4.x known constraints:**
- Pixel coordinates in `-at` mode must be **integers** — decimals cause `invalid arguments`
- `-autoring` is incompatible with `-at` mode (Siril bug — use fixed `setphot` values instead)

---

## 6. Selecting comparison stars

Comparison stars are fetched from **VizieR APASS DR9** (`II/336/apass9`) or via Siril's `findcompstars` command.

Automatic filters:
- `e_Vmag < 0.05` — catalog uncertainty < 0.05 mag
- `|Vmag − target_mag| < 2.0` — within 2 mag of target
- Separation from target > 5 arcsec — avoids blending

The pipeline selects up to 18 closest-in-magnitude stars.

---

## 7. Apparent magnitude computation

The pipeline computes an approximate apparent V magnitude using the APASS catalog magnitudes of the comparison stars:

```
V_app = V_C + ensemble_V
```

where:
- `V_C` is the differential magnitude from Siril's `light_curve` output
- `ensemble_V` is the **median** V catalog magnitude of the comparison stars (read from `comp_stars.csv`)

This is mathematically equivalent to Siril's AAVSO calibration formula:
```
V_std = (V_ins_target − V_ins_comp) + V_catalog_comp
```

**Limitations:** Because an LP filter is used instead of a standard V filter, `V_app` carries an additional color-dependent systematic offset of ±0.2–0.5 mag vs. the true V magnitude. Relative variations within a session are unaffected.

Output in `photometry.csv`:
```csv
# Ensemble V (APASS comp stars median): 12.284
# V_app = V_C + ensemble_V
JD,V_C,V_app,err
```

---

## 8. FWHM extraction

Per-frame seeing FWHM is extracted from Siril's registration sequence file (`.seq`).

The `.seq` file contains `R0` lines with registration star FWHM in pixels (x and y separately). These are correlated with the FITS frame timestamps (`DATE-OBS` header) and scaled by the plate scale (1.035 arcsec/px for the Seestar S30 Pro).

Output in `fwhm.csv`:
```csv
# FWHM per frame from Siril registration star detection
# Plate scale: 1.0350 arcsec/px
JD,FWHM_x_arcsec,FWHM_y_arcsec
```

Note: this is the registration star's FWHM (seeing proxy), not a PSF fit on the target star.

---

## 9. AAVSO Extended Format — field reference

```csv
#TYPE=EXTENDED
#OBSCODE=XXXX
#SOFTWARE=Siril+seestar_varstar_siril.py
#FILTER=CV
#DELIM=,
#DATE=JD
#OBSTYPE=CCD
NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
```

| Field | Value | Rationale |
|-------|-------|-----------|
| `FILT` | `CV` | LP filter is not a standard band; CV is correct for broadband/unfiltered |
| `TRANS` | `NO` | No color transformation applied |
| `MTYPE` | `DIFF` | Differential — relative to comparison ensemble |
| `CNAME` | `ENSEMBLE` | Multiple comparison stars |
| `CMAG` | `na` | Ensemble magnitude not reported as a single value |
| `AMASS` | `na` | Airmass not computed (requires observer coordinates) |
| `NOTES` | `LP_filter_Seestar_S30Pro` | Documents the filter |

> **Before submitting:** replace `OBSCODE=XXXX` with your real code from [aavso.org/register](https://www.aavso.org/register).

---

## 10. Post-processing: `varstar_postprod/`

After the pipeline, use these tools to produce publication-quality figures — independently of Siril.

### `plot_lightcurve.py`

```bash
python plot_lightcurve.py /path/to/results/StarName --bin 5 --sigma 3 --format both
```

| Flag | Default | Description |
|------|---------|-------------|
| `--bin` | 10 | Bin width in minutes |
| `--sigma` | 3.0 | Sigma-clipping threshold (MAD-based) |
| `--format` | both | `pdf`, `png`, or `both` |
| `--periodogram` | off | Add Lomb-Scargle periodogram |
| `--period` | — | Period (days) for phase-folded plot |
| `--t0` | — | Reference epoch (JD) for phase fold |
| `--no-diag` | off | Skip diagnostics figure |
| `--out` | same as input | Output directory |

Output files:
- `light_curve_pub.pdf/png` — light curve with inverse-variance weighted bins
- `diagnostics_pub.pdf/png` — error distribution and scatter quality
- `periodogram_pub.pdf/png` — Lomb-Scargle period search
- `phase_folded_pub.pdf/png` — if `--period` is given

### `analysis.ipynb`

Open with JupyterLab or VS Code. Edit the configuration cell:

```python
DATA_DIR = "/path/to/results/StarName"
SIGMA    = 3.0
BIN_MIN  = 10
PERIOD   = None   # set to period in days to enable phase folding
```

Sections: light curve · diagnostics · summary statistics · Lomb-Scargle · phase fold · FWHM vs time · LaTeX table row.

---

## 11. Siril commands reference (1.4 syntax)

```
link basename -out=directory
    Create a numbered sequence from all FITS in current directory.
    NOTE: -out= value must NOT be quoted in .ssf scripts.

register sequence -2pass
    Star-pattern matching registration. Offline, robust.

seqapplyreg sequence -framing=max -filter-round=Nk
    Apply registration transforms. -framing=max uses maximum field.

calibrate sequence [-bias=file] [-dark=file] [-flat=file] [-cc=banding]
    Apply calibration frames. -cc=banding corrects column banding (CMOS).

stack sequence [type] [lo] [hi] [-nonorm | -norm=mul] [-out=name]
    rej 3 3 = Winsorized sigma clipping 3σ. -nonorm for darks, -norm=mul for flats.

seqplatesolve sequence -nocache -force -focal=mm -pixelsize=um -radius=deg
    Per-frame plate solve using Gaia DR3 (online).

setphot -aperture=px -inner=px -outer=px -gain=value
    Configure aperture photometry. gain in e-/ADU.

findcompstars target_name -catalog=APASS -maxmag=N -out=comp_stars.csv
    Fetch APASS comparison stars around target.

light_curve sequence channel [-ninastars=csv] [-at=x,y] [-refat=x,y ...]
    Aperture photometry on every frame.
    channel: 0 for 1-layer mono.
    -ninastars: use a comp star CSV from findcompstars.
```

Full reference: [siril.readthedocs.io/en/latest/Commands.html](https://siril.readthedocs.io/en/latest/Commands.html)

---

## 12. Typical precision and observing strategy

### Expected precision with Seestar S30 Pro (30 mm aperture)

| Target magnitude | Expected σ per frame | Notes |
|-----------------|---------------------|-------|
| < 10 mag | > 0.05 mag | Near saturation — shorten exposure |
| 10–13 mag | 0.01–0.04 mag | Good regime |
| 13–14 mag | 0.05–0.10 mag | Usable for large-amplitude variables |
| > 14 mag | > 0.10 mag | Marginal — add calibration frames |

Scintillation sets a hard floor of ~20–40 mmag on a 30 mm aperture regardless of target brightness. Exoplanet transit detection (< 10 mmag) is not feasible with this instrument.

### Short-period variables (RR Lyrae, δ Scuti, eclipsing binaries)

- Use raw 20s frames — don't stack
- Cover at least 1–2 full periods
- Cadence: one point every ~25s (frame + readout)

### Long-period variables (Mira, SR, semi-regular)

- Stack groups of 5–10 frames to reduce scatter
- 1–3 measurements per night is sufficient
- Use **"Photometry only"** on repeat nights

### Observing tips for best precision

- **Observe near meridian** — lowest airmass, minimum extinction and scintillation
- **Good transparency** — avoid nights with haze or high humidity
- **No LP filter** (if dark sky allows) — removes color systematics

---

## 13. Troubleshooting

Run `test_connections.py` from Siril first — it checks all APIs and the Siril connection.

| Error | Cause | Fix |
|-------|-------|-----|
| `No sequence light_ found` | Path has spaces | Move session to a path without spaces |
| `PSF cannot be computed on channel 1` | Sequence is 1-layer mono | Fixed: script uses channel 0 |
| `plate solve had errors` | Too few stars or no internet | Ensure internet; try longer exposures |
| `light_curve.dat not found` | Target outside field or no WCS | Check coordinates; verify plate solve |
| `VSX: no results` | VizieR unreachable | Run `test_connections.py` |
| `invalid arguments` in light_curve | Decimal pixel coordinates | Fixed: pipeline rounds to int |

---

## References

- [Siril documentation](https://siril.readthedocs.io/en/latest/)
- [Siril photometry tutorial](https://siril.org/tutorials/photometry/)
- [AAVSO Extended File Format](https://www.aavso.org/aavso-extended-file-format)
- [AAVSO CCD photometry guide](https://www.aavso.org/ccd-photometry-guide)
- [AAVSO Light Curve Generator](https://app.aavso.org/lcg/)
- [VizieR VSX catalog](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=B/vsx/vsx)
- [VizieR APASS DR9](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/336/apass9)
