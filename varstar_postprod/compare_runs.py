#!/usr/bin/env python3
"""
compare_runs.py — Compare pipeline runs of the same star (e.g. two comparison-star sets).

Usage
-----
    python compare_runs.py <run_dir> <run_dir> [<run_dir> …] [options]

Each <run_dir> is a results/<star>/ directory (photometry.csv, pipeline.log).
The first run is the reference of the difference panel.

Options
-------
    --labels L …   Legend labels (default: the folder names)
    --format FMT   pdf, png or both (default: both)
    --out DIR      Output directory (default: the parent of the first run)
    --stem NAME    Output file name (default: compare_runs)

Figure
------
    top     V (or V − C) against time for every run, with error bars
    bottom  run − reference on the frames both runs measured, with a linear fit

Printed per run: points, comparison stars and their mean catalogue V (from
pipeline.log), point-to-point scatter (robust σ of successive differences / √2)
and median error bar; per difference: median offset, robust σ and the drift of
the linear fit over the night.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_lightcurve import RCPARAMS, load_photometry, save_fig  # noqa: E402

COLORS = ["#2471a3", "#c0392b", "#27ae60", "#8e44ad", "#d35400"]
ERR_COLORS = ["#9ec3e0", "#e8a9a1", "#a9dfbf", "#d2b4de", "#f5cba7"]
_ENSEMBLE = re.compile(r"Ensemble: (\d+) comparison stars?, flux-mean catalogue V = ([\d.]+)")


def robust_sigma(x: np.ndarray) -> float:
    """1.4826 × median absolute deviation."""
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def ensemble_info(run_dir: Path) -> tuple[int | None, float | None]:
    """Number of comparison stars and their mean catalogue V, from pipeline.log."""
    log = run_dir / "pipeline.log"
    if log.exists():
        for m in _ENSEMBLE.finditer(log.read_text(encoding="utf-8", errors="replace")):
            return int(m.group(1)), float(m.group(2))
    return None, None


def load_run(run_dir: Path, label: str) -> dict:
    d = load_photometry(run_dir / "photometry.csv")
    y = d["vapp"] if d["has_vapp"] else d["vc"]
    n_comp, comp_v = ensemble_info(run_dir)
    return {
        "label": label, "jd": d["jd"], "y": y, "err": d["err"],
        "has_vapp": d["has_vapp"], "star": d["star_name"],
        "n_comp": n_comp, "comp_v": comp_v,
        "p2p": robust_sigma(np.diff(y)) / np.sqrt(2) if y.size > 2 else float("nan"),
        "med_err": float(np.median(d["err"])),
    }


def difference(ref: dict, run: dict) -> dict | None:
    """run − ref on the frames both measured (matched on JD to 1e-5 d ≈ 1 s)."""
    _, ia, ib = np.intersect1d(np.round(ref["jd"], 5), np.round(run["jd"], 5),
                               return_indices=True)
    if ia.size < 3:
        return None
    jd = ref["jd"][ia]
    diff = run["y"][ib] - ref["y"][ia]
    x = jd - jd.min()
    slope, icpt = np.polyfit(x, diff, 1)
    return {"jd": jd, "diff": diff, "fit": slope * x + icpt,
            "median": float(np.median(diff)), "sigma": robust_sigma(diff),
            "drift": float(slope * np.ptp(x)), "hours": float(np.ptp(x) * 24), "n": ia.size}


def run_label(r: dict) -> str:
    comps = (f"{r['n_comp']} comps, ⟨V⟩ = {r['comp_v']:.2f}" if r["n_comp"]
             else "comps: n/a")
    return f"{r['label']}  ({comps};  σ$_{{p2p}}$ = {r['p2p']:.3f})"


def plot_runs(runs: list[dict], diffs: list[dict | None]) -> plt.Figure:
    quantity = "V" if all(r["has_vapp"] for r in runs) else "V − C"
    off = int(min(r["jd"].min() for r in runs))
    with plt.rc_context(RCPARAMS):
        fig, (ax, ad) = plt.subplots(
            2, 1, figsize=(7.6, 6.2), sharex=True,
            gridspec_kw={"height_ratios": [3, 1.4], "hspace": 0.06})
        for i, r in enumerate(runs):
            c, ce = COLORS[i % len(COLORS)], ERR_COLORS[i % len(ERR_COLORS)]
            ax.errorbar(r["jd"] - off, r["y"], yerr=r["err"], fmt="o", ms=3.2,
                        color=c, ecolor=ce, elinewidth=0.6, capsize=0,
                        alpha=0.9 if i else 1.0, zorder=3 + i, label=run_label(r))
        ax.invert_yaxis()
        ax.set_ylabel(f"{quantity}  [mag]")
        ax.legend(loc="best", fontsize=8.5)
        ax.set_title(f"{runs[0]['star']} — pipeline runs compared", fontweight="bold")

        ad.axhline(0.0, color="0.5", lw=0.8, zorder=1)
        for i, (r, d) in enumerate(zip(runs[1:], diffs), start=1):
            if d is None:
                continue
            c = COLORS[i % len(COLORS)]
            ad.plot(d["jd"] - off, d["diff"], "o", ms=2.6, color=c, alpha=0.8, zorder=3,
                    label=f"{r['label']} − {runs[0]['label']}:  median {d['median']:+.3f}, "
                          f"σ = {d['sigma']:.3f}, drift {d['drift']:+.3f} mag")
            ad.plot(d["jd"] - off, d["fit"], color=c, lw=1.2, zorder=4)
        ad.invert_yaxis()
        ad.set_ylabel("Δ  [mag]")
        ad.set_xlabel(f"JD − {off:,}  [d]")
        ad.legend(loc="best", fontsize=8)
    return fig


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Compare pipeline runs of the same star.")
    ap.add_argument("runs", nargs="+", type=Path, help="results/<star>/ directories")
    ap.add_argument("--labels", nargs="+", default=None, help="legend labels")
    ap.add_argument("--format", "-f", dest="fmt", choices=["pdf", "png", "both"],
                    default="both")
    ap.add_argument("--out", "-o", type=Path, default=None)
    ap.add_argument("--stem", default="compare_runs")
    args = ap.parse_args(argv)

    dirs = [p.expanduser().resolve() for p in args.runs]
    if len(dirs) < 2:
        print("[error] give at least two runs", file=sys.stderr)
        return 1
    for d in dirs:
        if not (d / "photometry.csv").exists():
            print(f"[error] no photometry.csv in {d}", file=sys.stderr)
            return 1
    labels = args.labels or [d.name for d in dirs]
    if len(labels) != len(dirs):
        print("[error] one label per run", file=sys.stderr)
        return 1

    runs  = [load_run(d, lab) for d, lab in zip(dirs, labels)]
    diffs = [difference(runs[0], r) for r in runs[1:]]

    print(f"{'run':28s} {'N':>4s} {'comps':>5s} {'<V>comp':>8s} {'σ_p2p':>7s} {'med err':>8s}")
    for r in runs:
        print(f"{r['label']:28s} {r['jd'].size:4d} {r['n_comp'] or '-':>5} "
              f"{r['comp_v'] if r['comp_v'] is not None else float('nan'):8.3f} "
              f"{r['p2p']:7.4f} {r['med_err']:8.4f}")
    for r, d in zip(runs[1:], diffs):
        if d is None:
            print(f"{r['label']} − {runs[0]['label']}: fewer than 3 common frames")
        else:
            print(f"{r['label']} − {runs[0]['label']}: {d['n']} common frames, "
                  f"median {d['median']:+.4f}, σ {d['sigma']:.4f}, "
                  f"drift {d['drift']:+.4f} mag over {d['hours']:.2f} h")

    out = (args.out or dirs[0].parent).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    for p in save_fig(plot_runs(runs, diffs), out, args.stem, args.fmt):
        print(f"saved: {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
