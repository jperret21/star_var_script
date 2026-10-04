#!/usr/bin/env python3
"""
rrlyrae_template_fit.py — Fit RR Lyrae light-curve templates to pipeline photometry.

The light curve is folded on the pulsation period and compared with the
empirical templates of Sesar et al. (2010, ApJ 708, 717), built from SDSS
Stripe 82 RR Lyrae. Each template T(φ) is normalised to 0 at maximum and 1 at
minimum light, with φ = 0 at maximum. For a fixed period P the model is

    m(t) = m_max + A · T( (t − T_max) / P  mod 1 )

with three free parameters: the magnitude at maximum m_max, the amplitude A
and the epoch of maximum T_max. Every template of the library is fitted and
the one with the lowest χ² is kept.

Usage
-----
    python rrlyrae_template_fit.py <path> [<path> …] [options]

<path> is a results/<star>/ directory or a photometry.csv. Several nights of
the same star can be given; they are fitted together.

Options
-------
    --period P      Pulsation period in days (default: VSX, via VizieR)
    --band B        Template band: u, g, r, i or z (default: r)
    --type T        Template library: ab, c or all (default: ab)
    --template ID   Fit this template only (e.g. 105)
    --amplitude A   Fix the amplitude (mag), e.g. from a night that covered
                    the whole rise — for a light curve covering only part of
                    the cycle
    --clip K        Reject points more than K σ from the fit (default: 4; 0: off)
    --bootstrap N   Block-bootstrap draws for the uncertainties (default: 200;
                    0: covariance matrix only)
    --radec RA DEC  Target coordinates in degrees (default: from pipeline.log)
    --templates P   Template archive or folder (default: download and cache)
    --format FMT    pdf, png or both (default: both)
    --out DIR       Output directory (default: the first input directory)
    --dark          Dark figure theme
    --no-diag       Skip the diagnostic figure

Outputs
-------
    rrlyrae_template_fit.{pdf,png}       folded light curve, template, residuals
    rrlyrae_template_fit_diag.{pdf,png}  Δχ² profile of T_max, χ²_ν per template
    rrlyrae_template_fit.json            fitted parameters and uncertainties

Notes
-----
* Times are converted from JD (UTC) to HJD (UTC) with a low-precision solar
  ephemeris (error of a few seconds, well below the timing precision of a
  template fit) when the target coordinates are known; otherwise the fit
  runs on JD.
* The templates are SDSS-band shapes. On V or unfiltered data the amplitude
  is free, but the shape differs slightly between bands; the residual panel
  shows it. g and r bracket V.
* Uncertainties come from a moving-block bootstrap that refits the whole
  library, so they include the choice of template and residuals correlated
  in time. The covariance-matrix errors (scaled by χ²_ν when it exceeds 1)
  and the Δχ² profile of T_max are reported too; they are smaller when the
  template does not match the data shape perfectly (χ²_ν > 1).
* With partial phase coverage, amplitude and T_max are correlated: the script
  reports the coverage and the A–T_max correlation, and flags T_max as
  extrapolated when no point lies near maximum light.
* O − C is given against the VSX ephemeris; check the time scale of the VSX
  epoch before publishing it.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import sys
import tarfile
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_lightcurve import RCPARAMS, load_photometry, save_fig  # noqa: E402

# ── Sesar et al. (2010) templates ─────────────────────────────────────────────
# The original URL (astro.washington.edu/users/ivezic/sdss/catalogs/) is gone;
# astroML mirrors the same archive.
TEMPLATE_URL    = ("https://github.com/astroML/astroML-data/raw/main/datasets/"
                   "RRLyr_ugriz_templates.tar.gz")
TEMPLATE_SHA256 = "f88cbdb35d6ba98945fae198c4e7c80aaadd22f0a132e10ad894e080838e7565"
TEMPLATE_CACHE  = Path.home() / ".cache" / "star_var_script" / "sesar2010"
TEMPLATE_REF    = "Sesar et al. (2010)"
_TEMPLATE_NAME  = re.compile(r"(\d+)([ugriz])\.dat")
RRC_IDS         = {0, 1}   # near-sinusoidal RRc shapes; 100+ are RRab

VIZIER       = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
LIGHT_DAY_AU = 0.0057755183          # light travel time over 1 au [d]

AMP_BOUNDS   = (0.05, 2.5)           # mag — keeps a fit to a short, nearly
                                     # linear segment from running away
N_GRID       = 2000                  # phase-shift grid (0.0005 in phase)
NEAR_MAX     = 0.03                  # |φ| below which maximum counts as observed
COVER_BINS   = 100                   # phase bins for the coverage fraction


# ═══════════════════════════════════════════════════════════════════════════════
# Templates
# ═══════════════════════════════════════════════════════════════════════════════

def _extract_templates(archive: bytes, dest: Path) -> int:
    """Write the <id><band>.dat members of the archive into dest; returns the count.

    Only regular files with a template name are written, under their base
    name, so the archive cannot write outside dest.
    """
    dest.mkdir(parents=True, exist_ok=True)
    n = 0
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            if member.isfile() and _TEMPLATE_NAME.fullmatch(name):
                (dest / name).write_bytes(tar.extractfile(member).read())
                n += 1
    return n


def template_dir(source: Path | None = None) -> Path:
    """Folder holding the template .dat files, downloading them on first use."""
    if source is not None and source.is_dir():
        return source
    if source is None and any(TEMPLATE_CACHE.glob("*.dat")):
        return TEMPLATE_CACHE
    if source is not None:
        archive = source.read_bytes()
    else:
        print(f"Downloading the {TEMPLATE_REF} templates …\n  {TEMPLATE_URL}")
        with urllib.request.urlopen(TEMPLATE_URL, timeout=60) as r:
            archive = r.read()
        digest = hashlib.sha256(archive).hexdigest()
        if digest != TEMPLATE_SHA256:
            raise RuntimeError(f"template archive checksum mismatch ({digest})")
    n = _extract_templates(archive, TEMPLATE_CACHE)
    if not n:
        raise RuntimeError("no template found in the archive")
    print(f"  {n} templates cached in {TEMPLATE_CACHE}")
    return TEMPLATE_CACHE


def load_templates(folder: Path, band: str, kind: str = "ab") -> dict[int, np.ndarray]:
    """Templates of one band as {id: array (N, 2) of phase, normalised mag}.

    kind: "ab", "c" or "all".
    """
    out: dict[int, np.ndarray] = {}
    for f in folder.glob(f"*{band}.dat"):
        m = _TEMPLATE_NAME.fullmatch(f.name)
        if not m or m.group(2) != band:
            continue
        tid = int(m.group(1))
        if kind == "all" or (kind == "c") == (tid in RRC_IDS):
            out[tid] = np.loadtxt(f)
    return dict(sorted(out.items()))


def template_type(tid: int) -> str:
    return "RRc" if tid in RRC_IDS else "RRab"


def template_mag(tpl: np.ndarray, phase: np.ndarray) -> np.ndarray:
    """Normalised template magnitude at any phase (periodic, linear interpolation)."""
    return np.interp(np.mod(phase, 1.0), tpl[:, 0], tpl[:, 1], period=1.0)


# ═══════════════════════════════════════════════════════════════════════════════
# Time and ephemeris
# ═══════════════════════════════════════════════════════════════════════════════

def hjd_correction(jd: np.ndarray, ra_deg: float, dec_deg: float) -> np.ndarray:
    """HJD − JD [d] for a star at (ra, dec), J2000.

    Low-precision solar coordinates (Astronomical Almanac): about 0.01° in
    longitude, i.e. a few seconds at most on the light travel time.
    """
    n   = np.asarray(jd, dtype=float) - 2451545.0
    L   = np.radians((280.460 + 0.9856474 * n) % 360.0)
    g   = np.radians((357.528 + 0.9856003 * n) % 360.0)
    lam = L + np.radians(1.915) * np.sin(g) + np.radians(0.020) * np.sin(2 * g)
    r   = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)
    eps = np.radians(23.439 - 0.0000004 * n)
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    # (unit vector to the star) · (geocentric Sun vector), equatorial frame
    dot = r * (np.cos(dec) * np.cos(ra) * np.cos(lam)
               + np.cos(dec) * np.sin(ra) * np.cos(eps) * np.sin(lam)
               + np.sin(dec) * np.sin(eps) * np.sin(lam))
    return -LIGHT_DAY_AU * dot


_TARGET_LINE = re.compile(r"^# Target\s*:.*RA=([-+\d.]+)\s+Dec=([-+\d.]+)")


def target_coords(star_dir: Path) -> tuple[float, float] | None:
    """Target (ra, dec) from the pipeline.log written next to photometry.csv."""
    log = star_dir / "pipeline.log"
    if not log.exists():
        return None
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        m = _TARGET_LINE.match(line)
        if m:
            return float(m.group(1)), float(m.group(2))
    return None


def vsx_ephemeris(ra: float, dec: float) -> dict | None:
    """Name, type, period and epoch of the nearest VSX star within 30″."""
    params = {"-source": "B/vsx/vsx", "-c": f"{ra:.6f} {dec:+.6f}", "-c.rs": "30",
              "-out": "Name,Type,Period,Epoch", "-out.add": "_r", "-sort": "_r"}
    try:
        with urllib.request.urlopen(f"{VIZIER}?{urllib.parse.urlencode(params)}",
                                    timeout=20) as r:
            text = r.read().decode("utf-8", errors="replace")
    except OSError:
        return None
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    try:
        header = [c.strip() for c in lines[0].split("\t")]
        first  = next(i for i, ln in enumerate(lines) if ln.startswith("-")) + 1
        row    = dict(zip(header, (c.strip() for c in lines[first].split("\t"))))
    except (IndexError, StopIteration):
        return None

    def num(key):
        try:
            return float(row.get(key, ""))
        except ValueError:
            return None

    return {"name": row.get("Name", ""), "type": row.get("Type", ""),
            "period": num("Period"), "epoch": num("Epoch")}


# ═══════════════════════════════════════════════════════════════════════════════
# Fit
# ═══════════════════════════════════════════════════════════════════════════════

def _grid_chi2(u: np.ndarray, y: np.ndarray, w: np.ndarray, tpl: np.ndarray,
               phi0: np.ndarray, amplitude: float | None
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """χ², m_max and A on a grid of phase shifts, (m_max, A) solved exactly.

    u: time in cycles, y: magnitudes (weighted mean removed), w: 1/σ².
    A is held inside AMP_BOUNDS (or fixed); m_max is then the exact optimum.
    """
    S, Sy, Syy = w.sum(), w @ y, w @ (y * y)
    chi2 = np.empty(phi0.size)
    m0   = np.empty(phi0.size)
    amp  = np.empty(phi0.size)
    for lo in range(0, phi0.size, 250):          # chunks keep N×G arrays small
        sl  = slice(lo, lo + 250)
        T   = template_mag(tpl, u[:, None] - phi0[None, sl])
        Sx  = w @ T
        Sxx = w @ (T * T)
        Sxy = (w * y) @ T
        if amplitude is None:
            D = S * Sxx - Sx ** 2
            with np.errstate(divide="ignore", invalid="ignore"):
                A = (S * Sxy - Sx * Sy) / D
            A = np.clip(np.nan_to_num(A, nan=AMP_BOUNDS[0]), *AMP_BOUNDS)
        else:
            A = np.full(Sx.shape, float(amplitude))
        m = (Sy - A * Sx) / S
        chi2[sl] = (Syy + m * m * S + A * A * Sxx
                    - 2 * m * Sy - 2 * A * Sxy + 2 * m * A * Sx)
        m0[sl], amp[sl] = m, A
    return chi2, m0, amp


def fit_template(u: np.ndarray, y: np.ndarray, err: np.ndarray, tpl: np.ndarray,
                 amplitude: float | None = None) -> dict:
    """Least-squares fit of one template to (u, y, err); u is time in cycles.

    Global search over a grid of phase shifts (m_max and A solved exactly at
    each), then a joint non-linear refinement of all parameters.
    Returns m_max, amplitude, phi0 (phase shift of maximum, in cycles), chi2,
    cov (over m_max, [A,] phi0, scaled by χ²_ν when > 1) and the grid profile.
    """
    from scipy.optimize import least_squares

    w    = 1.0 / err ** 2
    ybar = (w @ y) / w.sum()
    yc   = y - ybar
    grid = np.arange(N_GRID) / N_GRID
    chi2, m0, amp = _grid_chi2(u, yc, w, tpl, grid, amplitude)
    k = int(np.argmin(chi2))

    if amplitude is None:
        lo, hi = AMP_BOUNDS
        a0 = min(max(amp[k], lo + 1e-6), hi - 1e-6)
        x0 = [m0[k], a0, grid[k]]
        bounds = ([-np.inf, lo, grid[k] - 0.05], [np.inf, hi, grid[k] + 0.05])

        def resid(p):
            return (yc - p[0] - p[1] * template_mag(tpl, u - p[2])) / err
    else:
        x0 = [m0[k], grid[k]]
        bounds = ([-np.inf, grid[k] - 0.05], [np.inf, grid[k] + 0.05])

        def resid(p):
            return (yc - p[0] - amplitude * template_mag(tpl, u - p[1])) / err

    sol  = least_squares(resid, x0, bounds=bounds, x_scale="jac")
    p    = sol.x
    chi2_fit = float(sol.fun @ sol.fun)
    dof  = max(1, y.size - p.size)
    cov  = np.linalg.pinv(sol.jac.T @ sol.jac) * max(1.0, chi2_fit / dof)

    return {
        "m_max":     float(p[0] + ybar),
        "amplitude": float(p[1]) if amplitude is None else float(amplitude),
        "phi0":      float(p[-1]),
        "chi2":      chi2_fit,
        "dof":       dof,
        "cov":       cov,
        "grid":      grid,
        "grid_chi2": chi2,
    }


def fit_library(t: np.ndarray, y: np.ndarray, err: np.ndarray, period: float,
                templates: dict[int, np.ndarray], amplitude: float | None = None,
                clip: float = 4.0) -> dict:
    """Fit every template, with iterative σ-clipping of the best one's outliers.

    Returns the best fit (plus 'template', 'mask', 'all' — χ²_ν per template —
    and the reference time 't_ref' of the phase shift).
    """
    t_ref = float(np.floor(t.min()))
    u     = (t - t_ref) / period
    mask  = np.ones(t.size, dtype=bool)
    for _ in range(6):
        fits = {tid: fit_template(u[mask], y[mask], err[mask], tpl, amplitude)
                for tid, tpl in templates.items()}
        tid  = min(fits, key=lambda k: fits[k]["chi2"])
        best = fits[tid]
        if clip <= 0:
            break
        model = best["m_max"] + best["amplitude"] * template_mag(
            templates[tid], u - best["phi0"])
        scale = np.sqrt(max(1.0, best["chi2"] / best["dof"]))
        new   = np.abs(y - model) <= clip * err * scale
        if np.array_equal(new, mask) or new.sum() < 10:
            break
        mask = new
    best.update(template=tid, mask=mask, t_ref=t_ref,
                all={k: f["chi2"] / f["dof"] for k, f in fits.items()},
                all_chi2={k: f["chi2"] for k, f in fits.items()})
    return best


def _parabola_offset(c: np.ndarray, k: int) -> float:
    """Sub-grid offset of a minimum at index k, from a parabola through 3 points."""
    if 0 < k < c.size - 1:
        denom = c[k - 1] - 2 * c[k] + c[k + 1]
        if denom > 0:
            return float(np.clip(0.5 * (c[k - 1] - c[k + 1]) / denom, -0.5, 0.5))
    return 0.0


def bootstrap(t: np.ndarray, y: np.ndarray, err: np.ndarray, period: float,
              templates: dict[int, np.ndarray], fit: dict,
              amplitude: float | None = None, n_boot: int = 200,
              seed: int = 1) -> dict:
    """Moving-block bootstrap (Künsch 1989) of the whole library fit.

    Draws blocks of consecutive points (length ≈ N^(1/3)), so residuals
    correlated in time stay correlated, refits every template on a ±0.1-cycle
    window of phase shifts around the best fit and keeps the best template
    each time. The spread therefore includes the choice of template.
    Returns the draws: T_max offset from the best fit [d], m_max, A, template.
    """
    rng  = np.random.default_rng(seed)
    mask = fit["mask"]
    u    = (t[mask] - fit["t_ref"]) / period
    yy, ee = y[mask], err[mask]
    n    = u.size
    L    = max(1, round(n ** (1 / 3)))
    nb   = int(np.ceil(n / L))
    step = 1.0 / N_GRID
    window = fit["phi0"] + np.arange(-200, 201) * step
    draws = []
    for _ in range(n_boot):
        starts = rng.integers(0, n - L + 1, nb)
        idx = (starts[:, None] + np.arange(L)).ravel()[:n]
        ub, yb, eb = u[idx], yy[idx], ee[idx]
        w    = 1.0 / eb ** 2
        ybar = (w @ yb) / w.sum()
        best = None
        for tid, tpl in templates.items():
            c, m0, amp = _grid_chi2(ub, yb - ybar, w, tpl, window, amplitude)
            k = int(np.argmin(c))
            if best is None or c[k] < best[0]:
                phi = window[k] + _parabola_offset(c, k) * step
                best = (c[k], (phi - fit["phi0"]) * period, m0[k] + ybar, amp[k], tid)
        draws.append(best[1:])
    d_t, m_max, amp, tid = (np.array(x) for x in zip(*draws))
    return {"dt": d_t, "m_max": m_max, "amplitude": amp, "template": tid,
            "block": L, "n": n_boot}


def _spread(x: np.ndarray) -> float:
    """Half the 16–84 % range: a 1σ equivalent that tolerates skewed draws."""
    lo, hi = np.percentile(x, [15.865, 84.135])
    return float((hi - lo) / 2)


def coverage_segments(phase: np.ndarray, nbins: int = COVER_BINS
                      ) -> list[tuple[float, float]]:
    """Phase intervals [lo, hi) holding data, on a grid of nbins bins."""
    occ = np.zeros(nbins, dtype=bool)
    occ[np.floor(np.mod(phase, 1.0) * nbins).astype(int) % nbins] = True
    segs, start = [], None
    for i, o in enumerate(np.append(occ, False)):
        if o and start is None:
            start = i
        elif not o and start is not None:
            segs.append((start / nbins, i / nbins))
            start = None
    return segs


def _profile_interval(grid: np.ndarray, chi2: np.ndarray, scale: float,
                      phi_best: float) -> tuple[float, float, bool]:
    """Phase-shift interval where Δχ²/scale ≤ 1 around the best shift.

    Returns (lo, hi) relative to phi_best, in cycles, and whether the Δχ² ≤ 1
    region is split in several pieces (a multimodal fit).
    """
    d   = np.mod(grid - phi_best + 0.5, 1.0) - 0.5          # wrapped to [-0.5, 0.5)
    ok  = (chi2 - chi2.min()) / scale <= 1.0
    order = np.argsort(d)
    d, ok = d[order], ok[order]
    i0  = int(np.argmin(np.abs(d)))
    lo = hi = i0
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < d.size - 1 and ok[hi + 1]:
        hi += 1
    multimodal = bool(ok.sum() > hi - lo + 1)
    return float(d[lo]), float(d[hi]), multimodal


def summarise(fit: dict, t: np.ndarray, period: float, boot: dict | None = None) -> dict:
    """Epoch of maximum, uncertainties and fit-quality figures from a fit.

    The *_err values are the headline 1σ: from the bootstrap when it ran
    (see bootstrap()), else from the covariance matrix; both are kept.
    """
    mask  = fit["mask"]
    cov   = fit["cov"]
    free_a = cov.shape[0] == 3
    phase = np.mod((t - fit["t_ref"]) / period - fit["phi0"], 1.0)

    # Epoch of maximum: the one closest to the middle of the data
    e0     = fit["t_ref"] + fit["phi0"] * period
    t_mid  = float(np.median(t[mask]))
    epoch  = e0 + np.round((t_mid - e0) / period) * period
    s_phi  = float(np.sqrt(cov[-1, -1]))
    scale  = max(1.0, fit["chi2"] / fit["dof"])
    lo, hi, multi = _profile_interval(fit["grid"], fit["grid_chi2"], scale, fit["phi0"])

    near_max = np.minimum(phase[mask], 1.0 - phase[mask]) < NEAR_MAX
    cover    = len(np.unique(np.floor(phase[mask] * COVER_BINS)))
    ranked   = sorted(fit["all_chi2"], key=fit["all_chi2"].get)
    cov_err  = {"epoch_max": s_phi * period,
                "m_max": float(np.sqrt(cov[0, 0])),
                "amplitude": float(np.sqrt(cov[1, 1])) if free_a else None}
    if boot is not None:
        err = {"epoch_max": _spread(boot["dt"]), "m_max": _spread(boot["m_max"]),
               "amplitude": _spread(boot["amplitude"]) if free_a else None}
        ids, counts = np.unique(boot["template"], return_counts=True)
        boot_info = {"n": boot["n"], "block_length": boot["block"],
                     "template_frequency": {str(int(i)): int(c) / boot["n"]
                                            for i, c in zip(ids, counts)}}
    else:
        err, boot_info = cov_err, None
    return {
        "epoch_max":        float(epoch),
        "epoch_max_err":    err["epoch_max"],
        "epoch_max_dchi2_interval": [lo * period, hi * period],
        "epoch_max_multimodal": multi,
        "max_observed":     bool(near_max.any()),
        "m_max":            fit["m_max"],
        "m_max_err":        err["m_max"],
        "amplitude":        fit["amplitude"],
        "amplitude_err":    err["amplitude"],
        "amplitude_fixed":  not free_a,
        "error_method":     "block bootstrap" if boot is not None else "covariance",
        "covariance_err":   cov_err,
        "bootstrap":        boot_info,
        "corr_amplitude_epoch": (float(cov[1, 2] / np.sqrt(cov[1, 1] * cov[2, 2]))
                                 if free_a else None),
        "chi2":             fit["chi2"],
        "dof":              fit["dof"],
        "chi2_red":         fit["chi2"] / fit["dof"],
        "n_points":         int(mask.sum()),
        "n_clipped":        int((~mask).sum()),
        "phase_coverage":   cover / COVER_BINS,
        "template":         fit["template"],
        "template_type":    template_type(fit["template"]),
        "next_template":    ranked[1] if len(ranked) > 1 else None,
        "next_delta_chi2":  (fit["all_chi2"][ranked[1]] - fit["chi2"]
                             if len(ranked) > 1 else None),
        "chi2_red_per_template": {str(k): v for k, v in fit["all"].items()},
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Figures
# ═══════════════════════════════════════════════════════════════════════════════

DARK = {
    "figure.facecolor": "#0d0d0f", "axes.facecolor": "#0d0d0f",
    "savefig.facecolor": "#0d0d0f", "axes.edgecolor": "#7a7a80",
    "axes.labelcolor": "#e6e6e6", "text.color": "#e6e6e6",
    "xtick.color": "#b0b0b5", "ytick.color": "#b0b0b5",
    "axes.grid": True, "grid.color": "#26262b", "grid.linewidth": 0.6,
    "legend.facecolor": "#16161a", "legend.edgecolor": "#3a3a40",
}


def _style(dark: bool) -> dict:
    return {**RCPARAMS, **DARK} if dark else RCPARAMS


def plot_fit(data: dict, t: np.ndarray, y: np.ndarray, err: np.ndarray,
             fit: dict, res: dict, tpl: np.ndarray, period: float, band: str,
             time_label: str, quantity: str, dark: bool = False) -> plt.Figure:
    """Folded light curve over two cycles with the fitted template, and residuals."""
    mask   = fit["mask"]
    phase  = np.mod((t - fit["t_ref"]) / period - fit["phi0"], 1.0)
    model  = res["m_max"] + res["amplitude"] * template_mag(tpl, phase)
    t_off  = int(np.floor(t.min()))
    c_tpl  = "#d4a94a" if dark else "#1b1b1b"
    c_err  = "#5a5a60" if dark else "#b0b0b0"
    c_seg  = "#e8d9a8" if dark else "#9e9e9e"
    edge   = "#0d0d0f" if dark else "0.25"

    with plt.rc_context(_style(dark)):
        fig, (ax, axr) = plt.subplots(
            2, 1, figsize=(7.6, 5.8), sharex=True,
            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06})

        segs = coverage_segments(phase[mask])
        for a in (ax, axr):
            for lo, hi in segs:
                for k in (0, 1):
                    a.axvspan(lo + k, hi + k, color=c_seg, alpha=0.10 if dark else 0.15,
                              lw=0, zorder=0)

        ph = np.linspace(0.0, 2.0, 4001)
        ax.plot(ph, res["m_max"] + res["amplitude"] * template_mag(tpl, ph),
                color=c_tpl, lw=1.6, zorder=2)

        tcol = t[mask] - t_off
        for k in (0, 1):
            x = phase[mask] + k
            for a, v in ((ax, y[mask]), (axr, (y - model)[mask])):
                a.errorbar(x, v, yerr=err[mask], fmt="none", ecolor=c_err,
                           elinewidth=0.6, capsize=0, zorder=3)
                sc = a.scatter(x, v, c=tcol, cmap="viridis", vmin=tcol.min(),
                               vmax=max(tcol.max(), tcol.min() + 1e-6), s=16,
                               edgecolors=edge, linewidths=0.3, zorder=4)
            if (~mask).any():
                xo = phase[~mask] + k
                ax.scatter(xo, y[~mask], marker="x", s=18, color="#e74c3c",
                           linewidths=0.8, zorder=5)
                axr.scatter(xo, (y - model)[~mask], marker="x", s=18,
                            color="#e74c3c", linewidths=0.8, zorder=5)

        ax.invert_yaxis()
        ax.set_xlim(0.0, 2.0)
        ax.set_ylabel(f"{quantity}  [mag]")
        axr.axhline(0.0, color=c_tpl, lw=1.0, zorder=2)
        axr.invert_yaxis()
        axr.set_ylabel("O − C  [mag]")
        axr.set_xlabel(f"Phase   (P = {period:.6f} d,  0 = maximum light)")
        axr.xaxis.set_major_locator(ticker.MultipleLocator(0.25))
        r = np.abs((y - model)[mask]).max()
        axr.set_ylim(1.25 * r, -1.25 * r)

        cbar = fig.colorbar(sc, ax=[ax, axr], pad=0.02, fraction=0.04)
        cbar.set_label(f"{time_label} − {t_off:,}  [d]")

        handles = [
            plt.Line2D([], [], color=c_tpl, lw=1.6,
                       label=f"{TEMPLATE_REF} template {res['template']} "
                             f"({res['template_type']}, {band} band)"),
            Patch(facecolor=c_seg, alpha=0.3, label="Observed phase coverage"),
        ]
        if (~mask).any():
            handles.append(plt.Line2D([], [], marker="x", ls="", color="#e74c3c",
                                      label=f"Clipped ({res['n_clipped']})"))
        ax.legend(handles=handles, loc="best")

        e, se = res["epoch_max"], res["epoch_max_err"]
        ax.set_title(
            f"$T_{{max}}$ = {time_label} {e:.5f} ± {se:.5f}"
            f"{' (extrapolated)' if not res['max_observed'] else ''}",
            fontsize=9, pad=6)
        fig.suptitle(f"{data['star_name']} — RR Lyrae template fit",
                     fontweight="bold", y=0.985)
    return fig


def plot_diagnostics(fit: dict, res: dict, period: float, time_label: str,
                     dark: bool = False) -> plt.Figure:
    """Δχ² profile of the epoch of maximum, and χ²_ν of every template."""
    c_line = "#d4a94a" if dark else "#1b1b1b"
    c_ab   = "#3b7dd8"
    c_c    = "#27ae60"
    scale  = max(1.0, res["chi2_red"])
    with plt.rc_context(_style(dark)):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.0, 3.4),
                                     gridspec_kw={"width_ratios": [1.3, 1]})
        d = np.mod(fit["grid"] - fit["phi0"] + 0.5, 1.0) - 0.5
        o = np.argsort(d)
        dchi = (fit["grid_chi2"] - fit["grid_chi2"].min()) / scale
        a1.plot(d[o] * period * 1440.0, dchi[o], color=c_line, lw=1.2)
        a1.axhline(1.0, color="#e74c3c", lw=0.8, ls="--", label="Δχ² = 1")
        lo, hi = res["epoch_max_dchi2_interval"]
        a1.axvspan(lo * 1440.0, hi * 1440.0, color="#e74c3c", alpha=0.12, lw=0)
        span = max(30.0, 4 * max(abs(lo), abs(hi)) * 1440.0)
        a1.set_xlim(-span, span)
        a1.set_ylim(0, max(10.0, min(dchi.max(), 60.0)))
        a1.set_xlabel(f"$T_{{max}}$ − best  [min]")
        a1.set_ylabel(r"$\Delta\chi^2$ / max(1, $\chi^2_\nu$)")
        a1.set_title(f"Epoch of maximum — template {res['template']}", fontsize=10)
        a1.legend(loc="upper center")

        ids  = list(fit["all"])
        vals = [fit["all"][k] for k in ids]
        cols = [c_c if k in RRC_IDS else c_ab for k in ids]
        bars = a2.bar(range(len(ids)), vals, color=cols, width=0.75)
        bars[ids.index(res["template"])].set_edgecolor("#e74c3c")
        bars[ids.index(res["template"])].set_linewidth(1.5)
        a2.set_xticks(range(len(ids)))
        a2.set_xticklabels([str(k) for k in ids], rotation=90, fontsize=7)
        a2.xaxis.set_minor_locator(ticker.NullLocator())
        a2.set_yscale("log")
        a2.set_ylabel(r"$\chi^2_\nu$")
        a2.set_title("Template library", fontsize=10)
        a2.legend(handles=[Patch(color=c_ab, label="RRab"), Patch(color=c_c, label="RRc")]
                  if any(k in RRC_IDS for k in ids) else [Patch(color=c_ab, label="RRab")],
                  loc="upper right")
        fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

def _csv_path(p: Path) -> Path:
    return p / "photometry.csv" if p.is_dir() else p


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=f"Fit {TEMPLATE_REF} RR Lyrae templates to pipeline photometry.")
    ap.add_argument("paths", nargs="+", type=Path,
                    help="results/<star>/ directories or photometry.csv files")
    ap.add_argument("--period", "-p", type=float, default=None,
                    help="pulsation period [d] (default: VSX)")
    ap.add_argument("--band", "-b", choices=list("ugriz"), default="r",
                    help="template band (default: r)")
    ap.add_argument("--type", dest="kind", choices=["ab", "c", "all"], default="ab",
                    help="template library (default: ab)")
    ap.add_argument("--template", type=int, default=None,
                    help="fit this template only")
    ap.add_argument("--amplitude", type=float, default=None,
                    help="fix the amplitude [mag]")
    ap.add_argument("--clip", type=float, default=4.0,
                    help="σ-clipping threshold on the residuals (0: off)")
    ap.add_argument("--bootstrap", type=int, default=200,
                    help="block-bootstrap draws for the uncertainties (0: off)")
    ap.add_argument("--seed", type=int, default=1, help="bootstrap random seed")
    ap.add_argument("--radec", type=float, nargs=2, metavar=("RA", "DEC"),
                    help="target coordinates [deg] (default: from pipeline.log)")
    ap.add_argument("--templates", type=Path, default=None,
                    help="template archive (.tar.gz) or folder of .dat files")
    ap.add_argument("--format", "-f", dest="fmt", choices=["pdf", "png", "both"],
                    default="both")
    ap.add_argument("--out", "-o", type=Path, default=None)
    ap.add_argument("--dark", action="store_true", help="dark figure theme")
    ap.add_argument("--no-diag", action="store_true", help="skip the diagnostic figure")
    args = ap.parse_args(argv)

    # ── Photometry ────────────────────────────────────────────────────────────
    csvs = [_csv_path(p.expanduser().resolve()) for p in args.paths]
    for c in csvs:
        if not c.exists():
            print(f"[error] {c} not found", file=sys.stderr)
            return 1
    sets = [load_photometry(c) for c in csvs]
    names = {d["star_name"] for d in sets}
    if len(names) > 1:
        print(f"[warning] several star names in the inputs: {', '.join(sorted(names))}")
    use_vapp = all(d["has_vapp"] for d in sets)
    quantity = "V" if use_vapp else "V − C"
    jd  = np.concatenate([d["jd"] for d in sets])
    y   = np.concatenate([d["vapp"] if use_vapp else d["vc"] for d in sets])
    err = np.concatenate([d["err"] for d in sets])
    ok  = np.isfinite(jd) & np.isfinite(y) & np.isfinite(err) & (err > 0)
    jd, y, err = jd[ok], y[ok], err[ok]
    order = np.argsort(jd)
    jd, y, err = jd[order], y[order], err[order]
    data = sets[0]
    if y.size < 10:
        print(f"[error] only {y.size} valid points", file=sys.stderr)
        return 1

    # ── Coordinates → HJD ─────────────────────────────────────────────────────
    radec = tuple(args.radec) if args.radec else target_coords(csvs[0].parent)
    if radec:
        t, time_label = jd + hjd_correction(jd, *radec), "HJD"
    else:
        t, time_label = jd, "JD"
        print("[warning] target coordinates unknown — fitting on JD, not HJD "
              "(give --radec RA DEC)")

    # ── Period ────────────────────────────────────────────────────────────────
    vsx = vsx_ephemeris(*radec) if radec else None
    if args.period:
        period, p_source = args.period, "user"
    elif vsx and vsx["period"]:
        period, p_source = vsx["period"], f"VSX ({vsx['name']}, {vsx['type']})"
    else:
        print("[error] no period: give --period (VSX lookup failed or no coordinates)",
              file=sys.stderr)
        return 1

    # ── Templates ─────────────────────────────────────────────────────────────
    try:
        folder = template_dir(args.templates)
    except (OSError, RuntimeError) as e:
        print(f"[error] templates: {e}", file=sys.stderr)
        return 1
    templates = load_templates(folder, args.band, "all" if args.template is not None
                               else args.kind)
    if args.template is not None:
        if args.template not in templates:
            print(f"[error] no template {args.template} in band {args.band} "
                  f"(available: {', '.join(map(str, templates))})", file=sys.stderr)
            return 1
        templates = {args.template: templates[args.template]}
    if not templates:
        print(f"[error] no {args.kind} template in band {args.band}", file=sys.stderr)
        return 1

    # ── Fit ───────────────────────────────────────────────────────────────────
    fit  = fit_library(t, y, err, period, templates, args.amplitude, args.clip)
    boot = (bootstrap(t, y, err, period, templates, fit, args.amplitude,
                      args.bootstrap, args.seed) if args.bootstrap > 0 else None)
    res  = summarise(fit, t, period, boot)
    tpl = templates[res["template"]]
    nights = int(1 + (np.diff(t) > 0.3).sum())

    oc = None
    if vsx and vsx["period"] and vsx["epoch"]:
        n  = round((res["epoch_max"] - vsx["epoch"]) / vsx["period"])
        oc = {"vsx_period": vsx["period"], "vsx_epoch": vsx["epoch"], "cycle": n,
              "o_minus_c_d": res["epoch_max"] - (vsx["epoch"] + n * vsx["period"])}

    se_min = res["epoch_max_err"] * 1440.0
    ce     = res["covariance_err"]
    lo, hi = (x * 1440.0 for x in res["epoch_max_dchi2_interval"])
    print(f"\n{data['star_name']} — {TEMPLATE_REF} template fit "
          f"({args.band} band, {'template ' + str(args.template) if args.template is not None else args.kind + ' library'})")
    print(f"  Data        : {res['n_points']} points ({res['n_clipped']} clipped"
          f"{f' at {args.clip:g}σ' if args.clip > 0 else ''}), {nights} night(s), "
          f"{quantity}, {time_label} (UTC)")
    print(f"  Period      : {period:.7f} d  [{p_source}]")
    print(f"  Coverage    : {res['phase_coverage']:.0%} of the cycle")
    nxt = (f"   next: {res['next_template']} (Δχ² = {res['next_delta_chi2']:.1f})"
           if res["next_template"] is not None else "")
    print(f"  Template    : {res['template']} ({res['template_type']})   "
          f"χ²_ν = {res['chi2_red']:.2f} ({res['chi2']:.1f} / {res['dof']}){nxt}")
    print(f"  Errors      : {res['error_method']}"
          + (f" ({boot['n']} draws, blocks of {boot['block']} points)" if boot else ""))
    print(f"  T_max       : {time_label} {res['epoch_max']:.5f} ± {res['epoch_max_err']:.5f}"
          f"  (± {se_min:.1f} min)"
          + ("" if res["max_observed"] else "  — extrapolated, no point near maximum"))
    print(f"                covariance ± {ce['epoch_max'] * 1440.0:.1f} min, "
          f"Δχ² ≤ 1: {lo:+.1f} / {hi:+.1f} min")
    print(f"  m_max       : {res['m_max']:.3f} ± {res['m_max_err']:.3f} mag"
          f"   (covariance ± {ce['m_max']:.3f})")
    if res["amplitude_fixed"]:
        print(f"  Amplitude   : {res['amplitude']:.3f} mag (fixed)")
    else:
        print(f"  Amplitude   : {res['amplitude']:.3f} ± {res['amplitude_err']:.3f} mag"
              f"   (covariance ± {ce['amplitude']:.3f}; "
              f"corr(A, T_max) = {res['corr_amplitude_epoch']:+.2f})")
    if boot:
        freq = sorted(res["bootstrap"]["template_frequency"].items(), key=lambda kv: -kv[1])
        print("  Bootstrap   : best template " + ", ".join(
            f"{k} {v:.0%}" for k, v in freq[:4]))
    print(f"  m_min       : {res['m_max'] + res['amplitude']:.3f} mag")
    if oc:
        print(f"  O − C (VSX) : {oc['o_minus_c_d'] * 1440.0:+.1f} min  "
              f"(cycle {oc['cycle']}, P = {oc['vsx_period']} d, E0 = {oc['vsx_epoch']})")
    if res["epoch_max_multimodal"]:
        print("  [warning] several separate T_max solutions within Δχ² ≤ 1 — "
              "see the diagnostic figure")
    if (not res["amplitude_fixed"] and res["phase_coverage"] < 0.5
            and abs(res["corr_amplitude_epoch"]) > 0.7):
        print("  [warning] partial phase coverage: amplitude and T_max are strongly "
              "correlated — fix the amplitude (--amplitude) if it is known")
    if not res["amplitude_fixed"] and np.isclose(res["amplitude"], AMP_BOUNDS, atol=1e-3).any():
        print(f"  [warning] amplitude at the fit bound {res['amplitude']} mag")

    # ── Outputs ───────────────────────────────────────────────────────────────
    out_dir = (args.out or csvs[0].parent).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = plot_fit(data, t, y, err, fit, res, tpl, period, args.band,
                   time_label, quantity, args.dark)
    saved = save_fig(fig, out_dir, "rrlyrae_template_fit", args.fmt)
    if not args.no_diag:
        fig = plot_diagnostics(fit, res, period, time_label, args.dark)
        saved += save_fig(fig, out_dir, "rrlyrae_template_fit_diag", args.fmt)

    report = {
        "star": data["star_name"], "inputs": [str(c) for c in csvs],
        "quantity": quantity, "time_scale": f"{time_label} (UTC)", "radec_deg": radec,
        "period_d": period, "period_source": p_source,
        "band": args.band, "library": args.kind if args.template is None else "single",
        "clip_sigma": args.clip, "nights": nights,
        "templates_source": {"reference": "Sesar et al. 2010, ApJ 708, 717",
                             "url": TEMPLATE_URL, "sha256": TEMPLATE_SHA256},
        **res, "o_minus_c_vsx": oc,
    }
    json_path = out_dir / "rrlyrae_template_fit.json"
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    saved.append(json_path)
    for p in saved:
        print(f"  saved: {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
