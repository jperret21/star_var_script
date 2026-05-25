# Pipeline — Technical Notes

Implementation decisions, known bugs worked around, and algorithm specifics. Not a tutorial.

---

## Version history

| Version | Key changes |
|---------|-------------|
| v0.1.5 | PRIMARY path via `findcompstars` + `-ninastars`; ensemble V_app; FWHM from `.seq` |
| v0.1.4 | Frame safety margin raised to 200 px |
| v0.1.3 | FALLBACK A/B comp star in-frame filtering |
| v0.1.2 | Early frame check in step 5; removed Siril annotate |
| v0.1.1 | WCS-based frame filtering |
| v0.1.0 | Initial release |

---

## v0.1.5 — Algorithm details

### Photometry path selection

The pipeline tries three strategies in order. The one that succeeds is logged as `[PRIMARY]`, `[FALLBACK A]`, or `[FALLBACK B]`.

**PRIMARY** — Siril `findcompstars` + `-ninastars`

```
findcompstars <target_name> -catalog=APASS -narrowband=0 -max_stars=N
  → comp_stars.csv   (type, name, ra, dec, mag)
light_curve r_light_ 0 -ninastars=comp_stars.csv
```

Siril handles star pixel lookup internally from the WCS headers. Requires WCS on every frame and at least 3 comp stars within the field.

**FALLBACK A** — pixel coordinates via VizieR query

```
light_curve r_light_ 0 -at=x,y -refat=x1,y1 -refat=x2,y2 ...
```

Comp stars queried from VizieR APASS DR9 (`e_Vmag < 0.05`, `|ΔV| < 2.0`), converted to pixel via TAN projection from the reference frame WCS. Coordinates are rounded to integers (Siril 1.4.x rejects decimals).

**FALLBACK B** — single nearest bright star

Last resort. One comp star only — differential magnitude is noisier and any intrinsic variability of the comp star is undetected.

---

### Apparent magnitude (V_app)

```
V_app = V_C + median(V_catalog[comp stars])
```

`V_C` is the raw Siril differential output (target − ensemble). `median(V_catalog)` is read from `comp_stars.csv` column 5 (Comp1 rows). This is mathematically the same as Siril's AAVSO formula but computed externally because `seqsetmag` is GUI-only and not scriptable.

**LP filter systematic:** both target and comp stars go through the same filter, so V_C is unaffected. The systematic enters through `median(V_catalog)` because APASS magnitudes are in standard V. For a CV (blue + red) vs. G/K comp stars, expect ±0.2–0.5 mag absolute offset. Relative variations within a session are reliable.

---

### FWHM extraction from `.seq`

Siril's registration sequence file has one `R0` line per frame containing the registration star FWHM in pixels (x, y separately). These are produced by the star-detection step of `register -2pass` — not a PSF fit on the target, but a good seeing proxy.

```
R0 <fwhm_x> <fwhm_y> <other fields...>
```

`R0` lines appear in the same order as the `I` (image) lines. The pipeline zips them 1:1, skips frames with `selected=0` or `fwhm_x ≤ 0`, reads `DATE-OBS` from each corresponding FITS header, converts to JD, and scales by plate scale (1.035 arcsec/px).

**Stem convention:** `registered.rstrip("_")` — e.g. `"r_light_"` → `"r_light"`. FITS filename is then `r_light_00012.fit` (stem + `_` + zero-padded index). If a frame is missing from disk (rejected by `seqapplyreg`), the I line still exists in the `.seq` with `selected=0`, so the zip alignment stays correct.

---

### Gain handling

`get_gain_eadu()` tries headers in order: `EGAIN`, `EPERDN`, `GAIN_E`, `CCDGAIN`, then `GAIN`. Values accepted only if `0.05 < g < 30`. The Seestar writes `GAIN=200` (ISO-equivalent, not e⁻/ADU), which is rejected by this range check.

**Current state:** always falls back to **1.0 e⁻/ADU**. The IMX585 at gain 200 is physically ~0.5–0.7 e⁻/ADU. Error bars are therefore ~20–40% wider than they should be, which is conservative and harmless for differential photometry.

**Fix when available:** the companion Seestar control app should write `EGAIN` with the real conversion factor into the FITS headers, and `get_gain_eadu()` will pick it up automatically.

---

## Siril 1.4.x constraints and workarounds

| Constraint | Workaround |
|-----------|------------|
| `-at` pixel coords must be integers | `round()` applied before building the command |
| `-autoring` incompatible with `-at` | Fixed `setphot` values used instead |
| `seqsetmag` is GUI-only, not scriptable | Ensemble V_app computed externally from `comp_stars.csv` |
| `.ssf` parser does not strip quotes from `-out="path"` | Always use `-out=path` (unquoted) |
| `seqpsf` in headless mode prints to console only, no file output | Not used; FWHM read from `.seq` instead |
| `light_curve -wcs` can fail if WCS availability check fails on some frames | PRIMARY path uses `-ninastars` (sky coords handled internally by Siril); FALLBACK A uses `-at` (bypasses WCS check) |

---

## Registration and alignment

`register -2pass` uses star-pattern matching (offline). The transform per frame is a **similarity** (translation + rotation + uniform scale) — no shear, no distortion correction.

`seqapplyreg -framing=max` pads frames to the union of all fields of view. This means the output frames are larger than the input. Stars near the border appear in fewer frames, which is why pixel coordinates outside the safety margin (200 px from each edge, v0.1.4+) are excluded from comp star selection.

`-filter-round=2.5k` keeps the 2500 best frames by star elongation. On a typical Seestar session of 200–500 frames this has no effect (all kept), but protects against including wind-shake or satellite trail frames.

---

## Plate scale

```
plate_scale = 206.265 × pixel_size_µm / focal_length_mm
            = 206.265 × 2.9 / 160
            ≈ 3.74 arcsec/px   (raw Bayer pixel)
```

After `seqapplyreg` the pixel scale is unchanged (Siril resamples to the same grid). The FWHM plate scale used for `.seq` → arcsec conversion is **1.035 arcsec/px** — this is the registered frame scale after the Seestar's internal preprocessing (the Seestar stacks 10s sub-exposures with slight drizzle, changing the effective pixel scale slightly).

> TODO: verify 1.035 empirically from a plate-solved frame (FITS CD matrix → arcsec/px).
