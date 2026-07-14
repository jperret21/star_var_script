# Seestar S30 Pro — Variable Star Photometry

Automated pipeline for variable star differential photometry with the **ZWO Seestar S30 Pro** (30 mm, f/5, Sony IMX585). Built on top of [Siril](https://siril.org) (open-source astronomical image processing) and driven entirely from its Script menu — no command line required.

---

## Features

- **Automatic target discovery** — queries VizieR/VSX to find all known variable stars in the field
- **Ensemble differential photometry** — APASS comparison stars fetched automatically from VizieR
- **Approximate apparent V magnitude** — computed from APASS catalog magnitudes of comp stars
- **Per-frame FWHM** — extracted from the Siril registration sequence and exported to CSV
- **AAVSO-ready export** — AAVSO Extended Format CSV, ready for WebObs submission
- **Publication-quality plots** — light curve, diagnostics, Lomb-Scargle periodogram, phase fold
- **Interactive notebook** — Jupyter notebook for in-depth session analysis

## Pipeline overview

```
Seestar FITS frames
    ↓
  Siril: calibrate → register → plate-solve → aperture photometry
    ↓
  photometry.csv  ·  fwhm.csv  ·  StarName_aavso.csv
    ↓
  varstar_postprod/plot_lightcurve.py  ·  analysis.ipynb

```

## Repository

```
star_var_script/
├── seestar_varstar_siril.py   Main script (run from Siril Script menu)
├── pipeline.py                Pipeline logic
├── varstar_postprod/          Post-processing tools (independent of Siril)
│   ├── plot_lightcurve.py     CLI plot generator
│   └── analysis.ipynb         Interactive analysis notebook
└── docs/
    └── photometry_guide.md    Full technical reference
```

## Quick start

1. Open Siril and set your session folder (the one containing `lights/`)
2. **Script → Run Script → `seestar_varstar_siril.py`**
3. Load the session, pick a variable star, fetch comp stars, click **Run**
4. Results are written to `results/StarName/`

After the pipeline completes, generate plots:

```bash
cd varstar_postprod
python plot_lightcurve.py /path/to/results/StarName --bin 5 --sigma 3 --format both
```

## Credits

Photometry pipeline powered by [Siril](https://siril.org) (© Siril team, GPLv3).
Catalog queries via [VizieR](https://vizier.cds.unistra.fr) (CDS Strasbourg) — VSX and APASS DR9.

## Documentation

Full technical reference (algorithms, step-by-step, output formats, Siril commands, AAVSO guide):
→ [docs/photometry_guide.md](docs/photometry_guide.md)
