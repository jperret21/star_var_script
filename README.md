# Seestar S30 Pro — Variable Star Photometry

Generate differential light curves from raw Seestar FITS frames, driven entirely from within Siril's Script menu. The pipeline discovers variable stars automatically via VizieR/VSX, aligns frames, plate-solves every frame for accurate sky coordinates, runs aperture photometry, computes approximate apparent V magnitudes, and exports an AAVSO-ready CSV.

---

## Repository layout

```
star_var_script/
├── seestar_varstar_siril.py     Main script (run from Siril Script menu)
├── pipeline.py                  Pure pipeline logic (no tkinter dependency)
├── test_connections.py          Connection diagnostics (run from Siril)
├── requirements.txt             Pinned deps for the .venv
├── docs/
│   └── photometry_guide.md      Technical reference
├── tools/                       Standalone CLI tools (need .venv)
│   ├── seestar_varstar_gui.py
│   └── seestar_variable_photometry.py
└── varstar_postprod/            Post-processing and analysis (independent of Siril)
    ├── plot_lightcurve.py       Publication-quality plot generator (CLI)
    └── analysis.ipynb           Interactive analysis notebook
```

---

## How it works

```
lights/*.fit
    │
    ├─ [optional] calibrate with darks / flats / bias
    │
    ├─ link + register -2pass        star-matching alignment (offline, robust)
    ├─ seqapplyreg                   apply transforms → r_light_*.fit
    ├─ seqplatesolve                 per-frame WCS from Gaia DR3
    │
    └─ light_curve (aperture phot)   differential magnitudes
         │
         ├─ photometry.csv           JD, V-C, V_app, err
         ├─ fwhm.csv                 JD, FWHM_x/y in arcsec (from .seq registration)
         └─ StarName_aavso.csv       AAVSO Extended Format
```

---

## Requirements

| Software | Version | Notes |
|----------|---------|-------|
| **Siril** | 1.4.3+ | [siril.org/download](https://siril.org/download/) — the script runs inside Siril |
| Internet | — | VizieR catalog queries (VSX + APASS) |

No extra Python packages needed — Siril ships its own Python 3.12 with `requests` included.

> For the virtual environment (`tools/` and `varstar_postprod/` scripts), see `requirements.txt`.

---

## Seestar S30 Pro hardware profile

| Parameter | Value | Source |
|-----------|-------|--------|
| Aperture | **30 mm** | |
| Focal length | 160 mm | FITS `FOCALLEN` |
| Pixel size | 2.9 µm | FITS `XPIXSZ` |
| Sensor | Sony IMX585 | GRBG Bayer pattern |
| Plate scale | ~1.035 arcsec/px | computed |
| Field of view | ~2.2° × 3.9° | computed |
| Gain (e⁻/ADU) | read from `EGAIN` header, fallback 1.0 | **not** the ISO value |
| Photometry channel | **0** (1-layer mono after registration) | — |

---

## Quick start

### 1. Open Siril and set your session directory

In Siril: **File → Open** or use the folder icon to navigate to your session folder (the one containing `lights/`).

### 2. Run the script

**Script → Run Script → `seestar_varstar_siril.py`**

A window opens. If it says "● Siril connected" in the top right, you're ready.

> If the script doesn't appear in the menu, make sure it is installed in:
> `~/Library/Application Support/org.siril.Siril/siril/scripts/`

### 3. Load the session

Type or browse to your session folder. Click **Load** — the script reads the FITS headers and shows the field center (RA / Dec) and number of frames.

### 4. Query variable stars

Click **↻ Query VSX catalog** — the table fills with all variable stars found in the field via VizieR/VSX.

### 5. Select a target

Click a star. Then click **Fetch comparison stars (APASS)** — the script queries VizieR/APASS for comparison stars.

### 6. Run

Click **▶ Prepare frames & Generate light curve**.

| Step | Operation | Duration (200 frames) |
|------|-----------|----------------------|
| 0 | Stack calibration masters (if any) | 1–3 min |
| 1+2 | link + register -2pass | 2–5 min |
| 3 | seqapplyreg | 3–8 min |
| 4 | seqplatesolve (Gaia DR3) | 15–40 min |
| 5 | setphot + light_curve | < 1 min |

---

## Resuming a partial run

Use the **Start from** dropdown:

| Option | When to use |
|--------|-------------|
| `Full pipeline (steps 1–5)` | Fresh session |
| `Plate solve + photometry (steps 4–5)` | Registration done |
| `Photometry only (step 5)` | Change target or comp stars and rerun instantly |

---

## Output files

```
session/
├── process/
│   ├── r_light_*.fit           Registered + plate-solved frames
│   ├── r_light_.seq            Siril registration sequence (used for FWHM)
│   ├── comp_stars.csv          APASS comparison stars (findcompstars output)
│   └── light_curve.dat         Raw Siril photometry output
└── results/
    └── StarName/
        ├── photometry.csv      JD, V-C, V_app, err
        ├── fwhm.csv            JD, FWHM_x_arcsec, FWHM_y_arcsec
        └── StarName_aavso.csv  AAVSO Extended Format
```

### `photometry.csv`

```csv
# Star: ES UMa
# Ensemble V (APASS comp stars median): 12.284
# V_app = V_C + ensemble_V  (approximate apparent V magnitude)
JD,V_C,V_app,err
2461161.334643,1.2421,13.5261,0.0868
...
```

`V_app` is an approximate apparent V magnitude computed as:
```
V_app = (V_target − V_comp_instrumental) + V_comp_catalog
```
This is mathematically equivalent to Siril's AAVSO formula. Because an LP filter is used instead of a standard V filter, `V_app` carries an additional color-dependent systematic offset of ±0.2–0.5 mag vs. the true V magnitude.

### `fwhm.csv`

Per-frame seeing FWHM extracted from the Siril registration sequence file (`.seq`), scaled by the plate scale:

```csv
# FWHM per frame from Siril registration star detection
# Plate scale: 1.0350 arcsec/px
JD,FWHM_x_arcsec,FWHM_y_arcsec
2461161.334643,4.520,4.310
...
```

### `StarName_aavso.csv`

[AAVSO Extended Format](https://www.aavso.org/aavso-extended-file-format), ready for WebObs submission. `FILT=CV` (LP filter), `MTYPE=DIFF`, `CNAME=ENSEMBLE`.

> **Before submitting:** replace `OBSCODE=XXXX` with your real observer code from [aavso.org](https://www.aavso.org/register).

---

## Post-processing: `varstar_postprod/`

After the pipeline runs, use these tools to produce publication-quality figures and interactive analysis — independently of Siril.

### `plot_lightcurve.py` — CLI plot generator

```bash
cd varstar_postprod
python plot_lightcurve.py /path/to/results/ES_UMa --bin 5 --sigma 3 --format both
```

Options:

| Flag | Default | Description |
|------|---------|-------------|
| `--bin` | 10 | Bin width in minutes |
| `--sigma` | 3.0 | Sigma-clipping threshold |
| `--format` | both | `pdf`, `png`, or `both` |
| `--periodogram` | off | Add Lomb-Scargle periodogram |
| `--period` | — | Period (days) for phase-folded plot |
| `--t0` | — | Reference epoch (JD) for phase fold |
| `--no-diag` | off | Skip diagnostics figure |
| `--out` | same as input | Output directory |

Output files (saved alongside `photometry.csv` or in `--out`):
- `light_curve_pub.pdf/png` — light curve with inverse-variance weighted bins
- `diagnostics_pub.pdf/png` — error distribution and scatter quality checks
- `periodogram_pub.pdf/png` — Lomb-Scargle period search
- `phase_folded_pub.pdf/png` — if `--period` is given

### `analysis.ipynb` — Interactive notebook

Open with JupyterLab or VS Code. Edit the configuration cell at the top:

```python
DATA_DIR = "/path/to/results/ES_UMa"
SIGMA    = 3.0
BIN_MIN  = 10
PERIOD   = None   # set to period in days to enable phase folding
```

Sections: light curve · diagnostics · summary statistics · Lomb-Scargle · phase fold · FWHM vs time · LaTeX table row.

---

## Scientific limitations

| Limitation | Impact | Mitigation |
|-----------|--------|-----------|
| LP filter ≠ standard V band | ±0.2–0.5 mag systematic on `V_app` | Use `FILT=CV` + note in AAVSO |
| No transformation applied | Color-dependent offset vs. standard system | Ensemble comp stars partially cancel |
| 30 mm aperture | Scintillation floor ~20–40 mmag; no exoplanet transit detection | Best suited for variables with amplitude > 0.2 mag |
| No calibration frames | Extra scatter from hot pixels, vignetting | Add darks + flats for sessions > 30 min |

---

## Troubleshooting

Run `test_connections.py` from Siril first — it checks all APIs and the Siril connection.

| Error | Cause | Fix |
|-------|-------|-----|
| `No sequence light_ found` | Path has spaces | Move session to a path without spaces |
| `PSF cannot be computed on channel 1` | Sequence is 1-layer mono | Fixed: script uses channel 0 |
| `plate solve had errors` | Too few stars or no internet | Ensure Siril has internet; try longer exposures |
| `light_curve.dat not found` | Target outside field or no WCS | Check coordinates; verify plate solve succeeded |
| `VSX: no results` | VizieR unreachable | Run `test_connections.py` |

---

## References

- [Siril documentation](https://siril.readthedocs.io/en/latest/)
- [Siril photometry tutorial](https://siril.org/tutorials/photometry/)
- [AAVSO Extended File Format](https://www.aavso.org/aavso-extended-file-format)
- [AAVSO CCD photometry guide](https://www.aavso.org/ccd-photometry-guide)
- [AAVSO Light Curve Generator](https://app.aavso.org/lcg/)
- [VizieR VSX catalog](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=B/vsx/vsx)
- [VizieR APASS DR9](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/336/apass9)
