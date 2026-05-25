# Pipeline — Technical Notes

---

## Version history

| Version | Changes |
|---------|---------|
| v0.1.5 | PRIMARY path via `findcompstars` + `-ninastars`; ensemble V_app; FWHM from `.seq` |
| v0.1.4 | Frame border safety margin raised to 200 px |
| v0.1.3 | In-frame filtering of comp stars for FALLBACK A/B |
| v0.1.2 | Early frame check at step 5; removed broken Siril annotate |
| v0.1.1 | WCS-based frame filtering |
| v0.1.0 | Initial release |

---

## v0.1.5 — How it works

### 1. Session path

The user provides the session folder path, which must contain a `lights/` subfolder with Seestar FITS frames. On load, the script reads the FITS headers of the first file to get:

- `OBJCTRA` / `OBJCTDEC` or `RA` / `DEC` → field center (equatorial coordinates)
- `NAXIS1` / `NAXIS2` + plate scale → field size in degrees
- Count of `.fit` files in `lights/`

### 2. VSX catalog query

A cone search is sent to **VizieR** (`B/vsx/vsx`) around the field center. Search radius is half the field diagonal, rounded up to the nearest 0.1° — wide enough to catch stars near the edges.

Default magnitude limit: **14**. All VSX variable stars within the cone and below this threshold are returned and displayed in the table. The user can filter by name or type in the UI and adjust the magnitude limit.

### 3. Target and comparison star selection

After selecting a target from the table, the script calls `findcompstars` (Siril) or queries VizieR APASS DR9 depending on availability. See [photometry path selection](#photometry-path-selection) for details.

### 4. Filter and calibration frames

**Filter:** choice between LP (Seestar default) or no filter. Only affects the `FILT` field in the AAVSO export (`CV` in both cases — neither matches a standard photometric band).

**Calibration:** the user can provide folders of darks, flats, and/or bias frames. If at least one type is given, the pipeline creates masters automatically at step 0:

```
link dark_ → stack dark_ rej 3 3 -nonorm  → master_dark.fit
link flat_  → stack flat_  rej 3 3 -norm=mul → master_flat.fit
link bias_  → stack bias_  rej 3 3 -nonorm  → master_bias.fit
calibrate light_ -dark=... -flat=... -bias=...  → pp_light_*.fit
```

Stacking uses Winsorized sigma clipping 3σ, robust against cosmic rays and satellite trails. Without calibration, the active sequence stays `light_`.

### 5. Run — steps and algorithms

**Steps 1+2 — Registration**

```
link light -out=process/
register light_ -2pass
```

`register -2pass` does offline star-pattern matching. The transform computed per frame is a similarity: translation + rotation + uniform scale — no distortion correction. Results (transform matrices) are written to `light_.seq`.

**Step 3 — Apply transforms**

```
seqapplyreg light_ -framing=max -filter-round=2.5k  →  r_light_*.fit
```

`-framing=max`: output frames use the union of all fields of view, so they are larger than the input. Stars near the edges appear in fewer frames — hence the 200 px border exclusion for comp star selection (v0.1.4).

`-filter-round=2.5k`: keeps the 2500 best frames by star elongation. On a typical Seestar session of 200–500 frames this usually rejects nothing, but guards against wind-shake or satellite trail frames.

**Step 4 — Per-frame plate solve**

```
seqplatesolve r_light_ -nocache -force -focal=160 -pixelsize=2.9 -radius=2.5
```

Each frame gets its own astrometric solution against the Gaia DR3 catalog (online). WCS headers (`CRVAL`, `CRPIX`, `CD` matrix) are written into every `r_light_*.fit`. Required because after `-framing=max`, the reference frame WCS does not apply to the others.

**Step 5 — Photometry**

```
findcompstars <target> -catalog=APASS -narrowband=0 -max_stars=N  →  comp_stars.csv
setphot -aperture=10 -inner=20 -outer=30 -gain=1.0
light_curve r_light_ 0 -ninastars=comp_stars.csv
```

`light_curve` runs aperture photometry on every frame. Channel `0` is the only layer in the mono frames produced by CFA registration. Raw output is `light_curve.dat` (JD, V-C, error).

**Python post-processing**

After Siril exits, the Python pipeline:
1. Reads `light_curve.dat`
2. Computes `V_app = V_C + median(V_catalog_comp)` from `comp_stars.csv`
3. Extracts per-frame FWHM from `r_light_.seq` (R0 lines), aligns with `DATE-OBS` from FITS headers, converts to arcsec
4. Writes result files

### 6. Output

```
results/StarName/
├── photometry.csv       JD, V_C, V_app, err
├── fwhm.csv             JD, FWHM_x_arcsec, FWHM_y_arcsec
└── StarName_aavso.csv   AAVSO Extended Format
```

`photometry.csv`:
```
# Ensemble V (APASS comp stars median): 12.284
# V_app = V_C + ensemble_V
JD,V_C,V_app,err
2461161.334643,1.2421,13.5261,0.0868
```

`fwhm.csv`:
```
# Plate scale: 1.0350 arcsec/px
JD,FWHM_x_arcsec,FWHM_y_arcsec
2461161.334643,4.520,4.310
```

---

## Photometry path selection

**PRIMARY** — `findcompstars` + `-ninastars`

Siril handles the sky-to-pixel coordinate conversion internally via the WCS headers. Requires valid WCS on all frames and at least 3 comp stars in the field.

**FALLBACK A** — pixel coordinates from VizieR

Comp stars queried from VizieR APASS DR9 (`e_Vmag < 0.05`, `|delta_V| < 2.0`), converted to pixels via TAN projection from the reference frame WCS. Coordinates rounded to integers — Siril 1.4.x rejects decimals.

**FALLBACK B** — single nearest bright star

Last resort. Differential magnitude is noisier, and any intrinsic variability of the comp star goes undetected.

---

## V_app and the LP filter bias

`V_app = V_C + median(V_catalog)`

V_C is unaffected by the LP filter (target and comp stars go through the same filter). The bias enters via `median(V_catalog)`, which is in standard APASS V while the observation is in LP. For a cataclysmic variable (blue + red) vs. G/K comp stars, expect ±0.2–0.5 mag absolute offset. Relative variations within a session are reliable.

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
206.265 * 2.9 µm / 160 mm = 3.74 arcsec/px   (raw Bayer pixel)
```

The value used for FWHM `.seq` → arcsec conversion is **1.035 arcsec/px** — the Seestar internally stacks 10s sub-exposures with slight drizzle, which shifts the effective pixel scale.

> TODO: verify 1.035 empirically from the CD matrix of a plate-solved frame.

---

## Gain e⁻/ADU

`get_gain_eadu()` tries headers in order: `EGAIN`, `EPERDN`, `GAIN_E`, `CCDGAIN`, then `GAIN`. Value accepted only if `0.05 < g < 30`. The Seestar writes `GAIN=200` (camera setting, not e⁻/ADU) — rejected by this range, falls back to **1.0 e⁻/ADU**.

Physical IMX585 value at gain 200: ~0.5–0.7 e⁻/ADU. Error bars are therefore ~20–40% too wide — conservative and harmless for differential photometry. The Seestar control app should write `EGAIN` into FITS headers; `get_gain_eadu()` will pick it up automatically.
