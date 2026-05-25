#!/usr/bin/env python3
"""
plot_lightcurve.py — Publication-quality light curves from Seestar variable star pipeline.

Usage
-----
    python plot_lightcurve.py <path> [options]

<path> can be:
    results/ES_UMa/      — directory of one star
    results/             — process all stars found inside
    /session/root/       — auto-discovers results/ subdirectory

Options
-------
    --period P     Period in days → enables phase-folded plot
    --t0 T0        Reference epoch (JD) for phase folding (default: first point)
    --bin N        Bin size in minutes (default: 5)
    --sigma N      Sigma-clipping threshold (default: 3.0)
    --format FMT   pdf, png, or both (default: both)
    --out DIR      Output directory (default: same as input data)
    --no-diag      Skip diagnostic figure
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec

# ── Publication-quality rcParams ──────────────────────────────────────────────
RCPARAMS: dict = {
    "font.family":          "serif",
    "font.size":            11,
    "axes.labelsize":       12,
    "axes.titlesize":       13,
    "xtick.labelsize":      10,
    "ytick.labelsize":      10,
    "xtick.direction":      "in",
    "ytick.direction":      "in",
    "xtick.top":            True,
    "ytick.right":          True,
    "xtick.minor.visible":  True,
    "ytick.minor.visible":  True,
    "xtick.major.size":     5.0,
    "xtick.minor.size":     2.5,
    "ytick.major.size":     5.0,
    "ytick.minor.size":     2.5,
    "axes.linewidth":       0.8,
    "legend.fontsize":      10,
    "legend.framealpha":    0.85,
    "legend.edgecolor":     "0.75",
    "figure.dpi":           150,
    "savefig.dpi":          300,
    "savefig.bbox":         "tight",
    "savefig.pad_inches":   0.05,
    "errorbar.capsize":     2.5,
    "lines.linewidth":      1.0,
}

# ── Colour palette ────────────────────────────────────────────────────────────
C_PTS     = "#2d2d2d"   # individual photometric points (dark grey)
C_ERR     = "#b0b0b0"   # error bar colour for individual pts
C_BIN_V   = "#c0392b"   # binned apparent magnitude (red)
C_BIN_VC  = "#2471a3"   # binned differential (blue)
C_CLIP    = "#e67e22"   # sigma-clipped outliers (orange)
C_PHASE   = "#7d3c98"   # phase-folded points (purple)


# ═══════════════════════════════════════════════════════════════════════════════
# Data loading
# ═══════════════════════════════════════════════════════════════════════════════

def load_photometry(path: Path) -> dict:
    """Load photometry.csv → dict with numpy arrays.

    Returns keys: jd, vc, vapp (None if not present), err, has_vapp,
                  ensemble_v, star_name, fwhm (None if not present), fwhm_jd.
    """
    jd, vc, vapp, err = [], [], [], []
    has_vapp  = False
    ensemble_v: float | None = None
    star_name:  str   | None = None

    with open(path, newline="", encoding="utf-8") as f:
        header_seen = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# Star:"):
                star_name = line.removeprefix("# Star:").strip()
            elif line.startswith("# Ensemble V"):
                try:
                    ensemble_v = float(line.split(":")[-1].strip())
                except ValueError:
                    pass
            elif line.startswith("#"):
                if "V_app" in line:
                    has_vapp = True
                continue
            elif not header_seen:
                header_seen = True
                if "V_app" in line:
                    has_vapp = True
                continue
            else:
                parts = line.split(",")
                try:
                    if has_vapp and len(parts) >= 4:
                        jd.append(float(parts[0]))
                        vc.append(float(parts[1]))
                        vapp.append(float(parts[2]))
                        err.append(float(parts[3]))
                    elif len(parts) >= 3:
                        jd.append(float(parts[0]))
                        vc.append(float(parts[1]))
                        err.append(float(parts[2]))
                except ValueError:
                    pass

    fwhm_path = path.parent / "fwhm.csv"
    fwhm_jd, fwhm_val = _load_fwhm(fwhm_path) if fwhm_path.exists() else (None, None)

    return {
        "jd":         np.array(jd,   dtype=float),
        "vc":         np.array(vc,   dtype=float),
        "vapp":       np.array(vapp, dtype=float) if vapp else None,
        "err":        np.array(err,  dtype=float),
        "has_vapp":   has_vapp and bool(vapp),
        "ensemble_v": ensemble_v,
        "star_name":  star_name or path.parent.name.replace("_", " "),
        "fwhm_jd":    fwhm_jd,
        "fwhm":       fwhm_val,
    }


def _load_fwhm(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load fwhm.csv → (jd_array, fwhm_arcsec_array)."""
    jd, fwhm = [], []
    with open(path, newline="", encoding="utf-8") as f:
        header_seen = False
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if not header_seen:
                header_seen = True
                continue
            parts = line.split(",")
            try:
                jd.append(float(parts[0]))
                fwhm.append(float(parts[1]))
            except (ValueError, IndexError):
                pass
    return np.array(jd, dtype=float), np.array(fwhm, dtype=float)


# ═══════════════════════════════════════════════════════════════════════════════
# Signal processing utilities
# ═══════════════════════════════════════════════════════════════════════════════

def sigma_clip(mag: np.ndarray, err: np.ndarray,
               nsigma: float = 3.0, maxiter: int = 10) -> np.ndarray:
    """Iterative sigma clipping on mag array. Returns boolean keep-mask."""
    mask = np.ones(len(mag), dtype=bool)
    for _ in range(maxiter):
        med = np.median(mag[mask])
        mad = np.median(np.abs(mag[mask] - med))
        std = 1.4826 * mad          # MAD → robust σ estimator
        if std == 0:
            break
        new_mask = np.abs(mag - med) <= nsigma * std
        if np.array_equal(new_mask, mask):
            break
        mask = new_mask
    return mask


def weighted_bin(jd: np.ndarray, mag: np.ndarray, err: np.ndarray,
                 bin_min: float = 5.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-variance weighted binning. Returns jd_bin, mag_bin, err_bin."""
    bin_days = bin_min / 1440.0
    edges    = np.arange(jd.min(), jd.max() + bin_days, bin_days)
    jd_b, mag_b, err_b = [], [], []
    for lo in edges[:-1]:
        sel = (jd >= lo) & (jd < lo + bin_days)
        if sel.sum() < 2:
            continue
        w   = 1.0 / err[sel] ** 2
        wt  = w.sum()
        jd_b.append(np.mean(jd[sel]))
        mag_b.append(np.dot(w, mag[sel]) / wt)
        err_b.append(1.0 / np.sqrt(wt))
    return np.array(jd_b), np.array(mag_b), np.array(err_b)


def phase_fold(jd: np.ndarray, period: float,
               t0: float | None = None) -> np.ndarray:
    """Return phase in [0, 1)."""
    if t0 is None:
        t0 = jd.min()
    return ((jd - t0) % period) / period


def lomb_scargle(jd: np.ndarray, mag: np.ndarray, err: np.ndarray,
                 pmin: float = 0.005, pmax: float = 2.0,
                 nfreq: int = 10_000) -> tuple[np.ndarray, np.ndarray, float]:
    """Lomb-Scargle periodogram via scipy. Returns (periods, power, best_period)."""
    from scipy.signal import lombscargle as _ls
    freqs  = np.linspace(1.0 / pmax, 1.0 / pmin, nfreq)
    ang    = 2 * np.pi * freqs
    t      = jd - jd.min()
    w      = 1.0 / err ** 2
    w     /= w.mean()
    power  = _ls(t, mag * w, ang, normalize=True)
    best   = 1.0 / freqs[np.argmax(power)]
    return 1.0 / freqs, power, best


# ═══════════════════════════════════════════════════════════════════════════════
# Formatting helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _jd_axis(ax: plt.Axes, jd: np.ndarray, label_bottom: bool = True) -> int:
    """Set x-axis limits and labels using JD − offset convention. Returns offset."""
    jd_offset = int(jd.min())
    ax.set_xlim(jd.min() - jd_offset - 0.002, jd.max() - jd_offset + 0.002)
    if label_bottom:
        ax.set_xlabel(f"JD−{jd_offset:,}  [d]", labelpad=4)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    return jd_offset


def _mag_axis(ax: plt.Axes, mag: np.ndarray, invert: bool = True,
              pad: float = 0.15) -> None:
    """Set y-axis limits for magnitude (inverted by default)."""
    lo, hi = mag.min() - pad, mag.max() + pad
    ax.set_ylim((hi, lo) if invert else (lo, hi))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))


def _stats_box(ax: plt.Axes, jd: np.ndarray, mag: np.ndarray, err: np.ndarray,
               n_clipped: int, sigma: float, loc: str = "lower left") -> None:
    """Add a statistics text box to an axes."""
    rms    = float(np.std(mag - np.median(mag)))
    merr   = float(np.median(err))
    txt    = (f"$N = {len(mag)}$\n"
              f"rms $= {rms:.3f}$ mag\n"
              f"$\\langle\\sigma\\rangle = {merr:.3f}$ mag")
    if n_clipped:
        txt += f"\n{n_clipped} rejected ({sigma}σ)"
    va  = "top"    if "upper" in loc else "bottom"
    xa  = 0.97     if "right" in loc else 0.03
    ya  = 0.97     if "upper" in loc else 0.03
    ha  = "right"  if "right" in loc else "left"
    ax.text(xa, ya, txt, transform=ax.transAxes, fontsize=9,
            verticalalignment=va, horizontalalignment=ha,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="0.75", alpha=0.92))


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 1 — main light curve
# ═══════════════════════════════════════════════════════════════════════════════

def plot_light_curve(data: dict, bin_min: float = 5.0,
                     sigma: float = 3.0) -> plt.Figure:
    """Main light curve: apparent V (top) + differential V-C (bottom) + optional FWHM."""
    with plt.rc_context(RCPARAMS):
        jd       = data["jd"]
        err      = data["err"]
        has_vapp = data["has_vapp"]
        star     = data["star_name"]
        has_fwhm = data["fwhm"] is not None

        n_panels = 2 if has_vapp else 1
        if has_fwhm:
            n_panels += 1

        heights = []
        if has_vapp:
            heights.append(3)          # apparent V
        heights.append(2)              # V-C (or only panel)
        if has_fwhm:
            heights.append(1.5)        # FWHM

        fig_h = 3.5 + 1.5 * n_panels
        fig   = plt.figure(figsize=(7.2, fig_h))
        gs    = GridSpec(n_panels, 1, hspace=0.08,
                         figure=fig, height_ratios=heights)

        panel_idx = 0
        ax_top    = None

        # ── Apparent magnitude panel ─────────────────────────────────────────
        if has_vapp:
            vapp = data["vapp"]
            mask = sigma_clip(vapp, err, sigma)
            n_cl = int((~mask).sum())

            ax = fig.add_subplot(gs[panel_idx])
            ax_top = ax
            panel_idx += 1

            off = int(jd.min())
            t   = jd - off

            # Clipped outliers (shown in a distinct colour for transparency)
            if n_cl:
                ax.errorbar(t[~mask], vapp[~mask], yerr=err[~mask],
                            fmt="x", color=C_CLIP, markersize=4,
                            elinewidth=0.6, ecolor=C_CLIP, alpha=0.55,
                            zorder=2, label=f"Rejected ({sigma}σ)")

            # Individual points
            ax.errorbar(t[mask], vapp[mask], yerr=err[mask],
                        fmt="o", color=C_PTS, markersize=2.5,
                        elinewidth=0.6, ecolor=C_ERR, alpha=0.55,
                        linewidth=0, zorder=3,
                        label=f"Individual ({mask.sum()} pts)")

            # Weighted bins
            t_b, v_b, e_b = weighted_bin(t[mask], vapp[mask], err[mask], bin_min)
            ax.errorbar(t_b, v_b, yerr=e_b,
                        fmt="s", color=C_BIN_V, markersize=6,
                        elinewidth=1.2, ecolor=C_BIN_V, linewidth=0,
                        zorder=5, label=f"{bin_min:.0f}-min bins")

            _mag_axis(ax, vapp[mask])
            ax.set_ylabel(r"$V$  [mag]", labelpad=4)
            _stats_box(ax, t[mask], vapp[mask], err[mask], n_cl, sigma,
                       loc="lower left")
            ax.legend(loc="upper right", markerscale=1.3)
            ax.set_title(star, fontsize=14, fontweight="bold", pad=6)
            plt.setp(ax.get_xticklabels(), visible=False)
            _jd_axis(ax, jd, label_bottom=False)

        # ── Differential V-C panel ───────────────────────────────────────────
        vc   = data["vc"]
        mask = sigma_clip(vc, err, sigma)
        n_cl = int((~mask).sum())

        share_x = ax_top if ax_top is not None else None
        ax_vc   = fig.add_subplot(gs[panel_idx], sharex=share_x)
        panel_idx += 1

        off = int(jd.min())
        t   = jd - off

        if n_cl:
            ax_vc.errorbar(t[~mask], vc[~mask], yerr=err[~mask],
                           fmt="x", color=C_CLIP, markersize=4,
                           elinewidth=0.6, ecolor=C_CLIP, alpha=0.55, zorder=2)

        ax_vc.errorbar(t[mask], vc[mask], yerr=err[mask],
                       fmt="o", color=C_PTS, markersize=2.5,
                       elinewidth=0.6, ecolor=C_ERR, alpha=0.55,
                       linewidth=0, zorder=3)

        t_b, v_b, e_b = weighted_bin(t[mask], vc[mask], err[mask], bin_min)
        ax_vc.errorbar(t_b, v_b, yerr=e_b,
                       fmt="s", color=C_BIN_VC, markersize=6,
                       elinewidth=1.2, ecolor=C_BIN_VC, linewidth=0, zorder=5)

        _mag_axis(ax_vc, vc[mask])
        ax_vc.set_ylabel(r"$V - C$  [mag]", labelpad=4)

        if ax_top is None:
            ax_top = ax_vc
            ax_vc.set_title(star, fontsize=14, fontweight="bold", pad=6)
            _stats_box(ax_vc, t[mask], vc[mask], err[mask], n_cl, sigma)
            ax_vc.legend(loc="upper right")

        # ── FWHM panel ───────────────────────────────────────────────────────
        if has_fwhm:
            fwhm_jd  = data["fwhm_jd"] - off
            fwhm_val = data["fwhm"]
            ax_fw    = fig.add_subplot(gs[panel_idx], sharex=ax_top)

            ax_fw.scatter(fwhm_jd, fwhm_val,
                          s=8, color="#27ae60", alpha=0.6,
                          linewidths=0, zorder=3)

            # Running median over 10 points
            if len(fwhm_val) >= 10:
                idx  = np.argsort(fwhm_jd)
                f_s  = fwhm_val[idx]
                kern = min(21, max(3, len(f_s) // 10) | 1)  # odd window
                from scipy.ndimage import median_filter
                f_sm = median_filter(f_s, size=kern)
                ax_fw.plot(fwhm_jd[idx], f_sm, color="#1a5e2e",
                           linewidth=1.2, zorder=4, label="Running median")

            ax_fw.set_ylabel(r"FWHM  [\"]", labelpad=4)
            ax_fw.set_ylim(bottom=0)
            ax_fw.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
            _jd_axis(ax_fw, jd, label_bottom=True)
            plt.setp(ax_vc.get_xticklabels(), visible=False)
        else:
            _jd_axis(ax_vc, jd, label_bottom=True)

        fig.align_ylabels()
        return fig


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 2 — data quality diagnostics
# ═══════════════════════════════════════════════════════════════════════════════

def plot_diagnostics(data: dict, sigma: float = 3.0) -> plt.Figure:
    """3-panel diagnostic figure: error histogram, mag–error scatter, mag histogram."""
    with plt.rc_context(RCPARAMS):
        jd       = data["jd"]
        err      = data["err"]
        has_vapp = data["has_vapp"]
        star     = data["star_name"]

        mag   = data["vapp"] if has_vapp else data["vc"]
        ylabel_lbl = r"$V$" if has_vapp else r"$V - C$"

        mask = sigma_clip(mag, err, sigma)
        mag_c = mag[mask]
        err_c = err[mask]

        fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.8))
        fig.suptitle(f"{star} — Photometric quality", fontsize=13,
                     fontweight="bold", y=1.01)

        # ── Error histogram ──────────────────────────────────────────────────
        ax = axes[0]
        ax.hist(err_c, bins=30, color=C_BIN_VC, edgecolor="white",
                linewidth=0.5, alpha=0.85)
        ax.axvline(np.median(err_c), color="#c0392b", linewidth=1.5,
                   linestyle="--", label=f"Median = {np.median(err_c):.3f} mag")
        ax.set_xlabel(r"Photometric error $\sigma$  [mag]", labelpad=4)
        ax.set_ylabel("Number of frames", labelpad=4)
        ax.legend(fontsize=9)
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

        # ── Magnitude vs error scatter ───────────────────────────────────────
        ax = axes[1]
        ax.scatter(err_c, mag_c, s=6, color=C_PTS, alpha=0.45,
                   linewidths=0, zorder=3)
        ax.invert_yaxis()
        ax.set_xlabel(r"Photometric error $\sigma$  [mag]", labelpad=4)
        ax.set_ylabel(f"{ylabel_lbl}  [mag]", labelpad=4)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))

        # ── Magnitude histogram ──────────────────────────────────────────────
        ax = axes[2]
        med  = np.median(mag_c)
        std  = np.std(mag_c)
        bins = np.linspace(mag_c.min() - 0.05, mag_c.max() + 0.05, 35)
        ax.hist(mag_c, bins=bins, color=C_BIN_V, edgecolor="white",
                linewidth=0.5, alpha=0.85, orientation="horizontal")
        ax.axhline(med, color="#2d2d2d", linewidth=1.5, linestyle="--",
                   label=f"Median = {med:.3f} mag")
        ax.axhline(med - std, color="#888", linewidth=1.0, linestyle=":",
                   label=f"±1σ = {std:.3f} mag")
        ax.axhline(med + std, color="#888", linewidth=1.0, linestyle=":")
        ax.invert_yaxis()
        ax.set_ylabel(f"{ylabel_lbl}  [mag]", labelpad=4)
        ax.set_xlabel("Number of frames", labelpad=4)
        ax.legend(fontsize=9)
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))

        fig.tight_layout()
        return fig


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 3 — phase-folded light curve
# ═══════════════════════════════════════════════════════════════════════════════

def plot_phase_folded(data: dict, period: float, t0: float | None = None,
                      bin_min: float = 1.0, sigma: float = 3.0) -> plt.Figure:
    """Phase-folded light curve. period in days."""
    with plt.rc_context(RCPARAMS):
        jd       = data["jd"]
        err      = data["err"]
        has_vapp = data["has_vapp"]
        star     = data["star_name"]

        mag    = data["vapp"] if has_vapp else data["vc"]
        ylabel = r"$V$  [mag]" if has_vapp else r"$V - C$  [mag]"

        mask  = sigma_clip(mag, err, sigma)
        phase = phase_fold(jd[mask], period, t0)
        mag_c = mag[mask]
        err_c = err[mask]

        # Extend by one cycle for visual clarity
        ph2   = np.concatenate([phase, phase + 1.0])
        m2    = np.concatenate([mag_c, mag_c])
        e2    = np.concatenate([err_c, err_c])

        # Bins in phase units
        bin_phase = bin_min / (period * 1440.0)
        ph_b, m_b, e_b = weighted_bin(ph2, m2, e2,
                                       bin_min=bin_phase * period * 1440.0)

        fig, ax = plt.subplots(figsize=(7.2, 4.5))

        ax.errorbar(ph2, m2, yerr=e2,
                    fmt="o", color=C_PTS, markersize=2.5, linewidth=0,
                    elinewidth=0.6, ecolor=C_ERR, alpha=0.45,
                    zorder=3, label=f"Individual ({mask.sum()} pts)")
        ax.errorbar(ph_b, m_b, yerr=e_b,
                    fmt="s", color=C_PHASE, markersize=7,
                    elinewidth=1.3, ecolor=C_PHASE, linewidth=0,
                    zorder=5, label=f"{bin_min:.0f}-min bins")

        ax.invert_yaxis()
        ax.set_xlabel(f"Phase  (P = {period:.4f} d)", labelpad=4)
        ax.set_ylabel(ylabel, labelpad=4)
        ax.set_xlim(-0.02, 2.02)
        ax.set_title(f"{star} — Phase-folded light curve", fontweight="bold")
        ax.axvline(1.0, color="0.6", linewidth=0.8, linestyle="--", alpha=0.7)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        ax.legend(loc="upper right")

        rms = float(np.std(mag_c - np.median(mag_c)))
        txt = (f"$P = {period:.4f}$ d = ${period * 1440:.1f}$ min\n"
               f"$N = {mask.sum()}$,  rms $= {rms:.3f}$ mag")
        if t0:
            txt = f"$T_0 = {t0:.4f}$\n" + txt
        ax.text(0.03, 0.96, txt, transform=ax.transAxes, fontsize=9,
                va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="0.75", alpha=0.92))
        return fig


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 4 — Lomb-Scargle periodogram
# ═══════════════════════════════════════════════════════════════════════════════

def plot_periodogram(data: dict, pmin: float = 0.005, pmax: float = 2.0,
                     sigma: float = 3.0) -> plt.Figure:
    """Lomb-Scargle periodogram (requires scipy)."""
    with plt.rc_context(RCPARAMS):
        jd       = data["jd"]
        err      = data["err"]
        has_vapp = data["has_vapp"]
        star     = data["star_name"]

        mag  = data["vapp"] if has_vapp else data["vc"]
        mask = sigma_clip(mag, err, sigma)

        periods, power, best = lomb_scargle(
            jd[mask], mag[mask], err[mask], pmin=pmin, pmax=pmax)

        fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.5),
                                 gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08})

        # Main periodogram
        ax = axes[0]
        ax.plot(periods, power, color="#2d2d2d", linewidth=0.7, alpha=0.85)
        ax.axvline(best, color=C_BIN_V, linewidth=1.4,
                   label=f"Best period: {best:.4f} d = {best * 1440:.1f} min")
        ax.set_ylabel("Normalised power", labelpad=4)
        ax.set_title(f"{star} — Lomb-Scargle periodogram",
                     fontweight="bold")
        ax.legend(loc="upper right")
        ax.set_xlim(pmin, pmax)
        ax.set_ylim(bottom=0)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        plt.setp(ax.get_xticklabels(), visible=False)

        # Log-scale insert for short periods
        ax2 = axes[1]
        ax2.plot(periods, power, color="#2d2d2d", linewidth=0.7, alpha=0.85)
        ax2.axvline(best, color=C_BIN_V, linewidth=1.4)
        ax2.set_xscale("log")
        ax2.set_xlabel("Period  [d]", labelpad=4)
        ax2.set_ylabel("Power", labelpad=4)
        ax2.set_xlim(pmin, pmax)
        ax2.set_ylim(bottom=0)
        ax2.xaxis.set_major_formatter(
            ticker.FuncFormatter(lambda x, _: f"{x:.3g}"))
        ax2.xaxis.set_minor_locator(ticker.LogLocator(subs=np.arange(2, 10) * 0.1,
                                                       numticks=12))
        ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))

        return fig


# ═══════════════════════════════════════════════════════════════════════════════
# Save helper
# ═══════════════════════════════════════════════════════════════════════════════

def save_fig(fig: plt.Figure, out_dir: Path, stem: str, fmt: str) -> list[Path]:
    """Save figure in requested format(s). Returns list of saved paths."""
    saved = []
    for f in (["pdf", "png"] if fmt == "both" else [fmt]):
        p = out_dir / f"{stem}.{f}"
        fig.savefig(p, format=f)
        saved.append(p)
    plt.close(fig)
    return saved


# ═══════════════════════════════════════════════════════════════════════════════
# Discovery
# ═══════════════════════════════════════════════════════════════════════════════

def find_star_dirs(root: Path) -> list[Path]:
    """Return list of directories containing photometry.csv."""
    if (root / "photometry.csv").exists():
        return [root]
    # results/ subdirectory
    if (root / "results").is_dir():
        root = root / "results"
    return sorted(p.parent for p in root.rglob("photometry.csv"))


# ═══════════════════════════════════════════════════════════════════════════════
# CLI entry point
# ═══════════════════════════════════════════════════════════════════════════════

def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description="Publication-quality light curves from Seestar pipeline data.")
    ap.add_argument("path", help="Results directory or session root")
    ap.add_argument("--period", "-p", type=float, default=None,
                    help="Orbital/pulsation period in days (enables phase plot + periodogram)")
    ap.add_argument("--t0", type=float, default=None,
                    help="Reference epoch JD for phase folding")
    ap.add_argument("--bin", "-b", dest="bin_min", type=float, default=5.0,
                    help="Bin size in minutes (default: 5)")
    ap.add_argument("--sigma", "-s", type=float, default=3.0,
                    help="Sigma-clip threshold (default: 3.0)")
    ap.add_argument("--format", "-f", dest="fmt",
                    choices=["pdf", "png", "both"], default="both",
                    help="Output format (default: both)")
    ap.add_argument("--out", "-o", type=Path, default=None,
                    help="Output directory (default: same as data)")
    ap.add_argument("--no-diag", action="store_true",
                    help="Skip diagnostic figure")
    ap.add_argument("--periodogram", action="store_true",
                    help="Compute and plot Lomb-Scargle periodogram")
    args = ap.parse_args(argv)

    root     = Path(args.path).expanduser().resolve()
    star_dirs = find_star_dirs(root)

    if not star_dirs:
        print(f"[error] No photometry.csv found under {root}", file=sys.stderr)
        sys.exit(1)

    for star_dir in star_dirs:
        csv_path = star_dir / "photometry.csv"
        out_dir  = args.out or star_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"[{star_dir.name}] Loading {csv_path.name} …")
        data = load_photometry(csv_path)

        if len(data["jd"]) == 0:
            print(f"  ⚠  No valid data points — skipping.")
            continue

        n  = len(data["jd"])
        vs = data["star_name"]
        ev = data["ensemble_v"]
        print(f"  Star       : {vs}")
        print(f"  Points     : {n}")
        print(f"  Ensemble V : {ev:.3f} mag" if ev else "  Ensemble V : n/a")
        if data["has_vapp"]:
            print(f"  V_app range: {data['vapp'].min():.3f} – {data['vapp'].max():.3f} mag")
        if data["fwhm"] is not None:
            print(f"  FWHM range : {data['fwhm'].min():.2f} – {data['fwhm'].max():.2f} \"")

        # Figure 1: light curve
        print("  → light curve …")
        fig = plot_light_curve(data, bin_min=args.bin_min, sigma=args.sigma)
        paths = save_fig(fig, out_dir, "light_curve_pub", args.fmt)
        for p in paths:
            print(f"     saved: {p}")

        # Figure 2: diagnostics
        if not args.no_diag:
            print("  → diagnostics …")
            fig = plot_diagnostics(data, sigma=args.sigma)
            paths = save_fig(fig, out_dir, "diagnostics_pub", args.fmt)
            for p in paths:
                print(f"     saved: {p}")

        # Figure 3: phase-folded
        if args.period:
            print(f"  → phase-folded (P={args.period} d) …")
            fig = plot_phase_folded(data, period=args.period, t0=args.t0,
                                    bin_min=max(1.0, args.bin_min / 5),
                                    sigma=args.sigma)
            paths = save_fig(fig, out_dir, "phase_folded_pub", args.fmt)
            for p in paths:
                print(f"     saved: {p}")

        # Figure 4: periodogram
        if args.periodogram or (args.period is None and n > 20):
            print("  → Lomb-Scargle periodogram …")
            try:
                fig = plot_periodogram(data, sigma=args.sigma)
                paths = save_fig(fig, out_dir, "periodogram_pub", args.fmt)
                for p in paths:
                    print(f"     saved: {p}")
            except ImportError:
                print("  ⚠  scipy not available — skipping periodogram")

        print()


if __name__ == "__main__":
    main()
