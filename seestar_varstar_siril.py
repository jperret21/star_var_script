#!/usr/bin/env python3
"""
Seestar S30 Pro — Variable Star Finder & Photometry
====================================================
Run from Siril: Script > Run Script > seestar_varstar_siril.py

Requirements:
  - Siril 1.4+ open, session folder set as working directory
  - Internet access for VizieR catalog queries

Troubleshoot:  Script > Run Script > test_connections.py
"""

from __future__ import annotations

import csv
import datetime
import math
import os
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
import tkinter as tk
from typing import Optional

# ── sirilpy — available when launched from Siril Script menu ─────────────────
try:
    sys.path.insert(0, "/Applications/Siril.app/Contents/Resources/share/siril/python_module")
    import sirilpy
    from sirilpy import SirilInterface
    HAS_SIRILPY = True
except ImportError:
    HAS_SIRILPY = False

# ── requests — bundled in Siril, also available in venv ──────────────────────
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# ─────────────────────────────────────────────────────────────────────────────
# Seestar S30 Pro hardware profile (from actual FITS headers)
# ─────────────────────────────────────────────────────────────────────────────
SEESTAR = {
    "focal":  160,    # mm  (FOCALLEN from FITS)
    "pixsz":  2.9,    # µm  (XPIXSZ from FITS)
    "fov_w":  2.2,    # deg width  (~2160px × 2.9µm / 160mm)
    "fov_h":  3.9,    # deg height (~3840px × 2.9µm / 160mm)
    # gain: read from FITS EGAIN header at runtime — GAIN=200 is ISO, not e-/ADU
}

VIZIER = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"

# AAVSO filter codes + note for each option presented to the user
# (code, extra_note)
FILTER_OPTIONS: dict[str, tuple[str, str]] = {
    "Clear — no filter (CV)":         ("CV",  ""),
    "LP anti-pollution filter (CV)":  ("CV",  "LP_filter_Seestar_S30Pro"),
    "V Johnson":                       ("V",   ""),
    "B Johnson":                       ("B",   ""),
    "R Johnson":                       ("R",   ""),
    "I Johnson":                       ("I",   ""),
    "Sloan r'":                        ("SRJ", ""),
    "Sloan i'":                        ("SIJ", ""),
}

# ─────────────────────────────────────────────────────────────────────────────
# Pure-Python FITS header reader (no astropy)
# ─────────────────────────────────────────────────────────────────────────────

def read_fits_header(path: Path) -> dict:
    """Read FITS keywords without astropy."""
    header: dict = {}
    try:
        with open(path, "rb") as f:
            while True:
                block = f.read(2880)
                if not block:
                    break
                for i in range(0, 2880, 80):
                    card = block[i:i + 80].decode("ascii", errors="replace")
                    if card.startswith("END"):
                        return header
                    if "=" in card[:10]:
                        key = card[:8].strip()
                        val = card[10:].split("/")[0].strip().strip("'").strip()
                        try:
                            header[key] = float(val)
                        except ValueError:
                            header[key] = val
    except Exception:
        pass
    return header


def sky_to_pixel(ra: float, dec: float, hdr: dict):
    """
    TAN gnomonic projection (ignores SIP, which is sub-pixel at Seestar scale):
    sky coordinates (J2000 degrees) → FITS pixel coordinates (1-based).
    Handles both CD-matrix and CDELT+PC formats.
    Returns (px, py) or (None, None) on failure.
    """
    try:
        crval1 = float(hdr["CRVAL1"])
        crval2 = float(hdr["CRVAL2"])
        crpix1 = float(hdr["CRPIX1"])
        crpix2 = float(hdr["CRPIX2"])
        if "CD1_1" in hdr:
            cd11 = float(hdr["CD1_1"])
            cd12 = float(hdr["CD1_2"])
            cd21 = float(hdr["CD2_1"])
            cd22 = float(hdr["CD2_2"])
        else:
            cdelt1 = float(hdr["CDELT1"])
            cdelt2 = float(hdr["CDELT2"])
            cd11 = cdelt1 * float(hdr.get("PC1_1", 1.0))
            cd12 = cdelt1 * float(hdr.get("PC1_2", 0.0))
            cd21 = cdelt2 * float(hdr.get("PC2_1", 0.0))
            cd22 = cdelt2 * float(hdr.get("PC2_2", 1.0))
        ra0   = math.radians(crval1)
        dec0  = math.radians(crval2)
        ra_r  = math.radians(ra)
        dec_r = math.radians(dec)
        dra   = ra_r - ra0
        denom = (math.sin(dec0) * math.sin(dec_r) +
                 math.cos(dec0) * math.cos(dec_r) * math.cos(dra))
        if abs(denom) < 1e-10:
            return None, None
        x = math.degrees(-math.cos(dec_r) * math.sin(dra) / denom)
        y = math.degrees((math.cos(dec0) * math.sin(dec_r) -
                          math.sin(dec0) * math.cos(dec_r) * math.cos(dra)) / denom)
        det = cd11 * cd22 - cd12 * cd21
        if abs(det) < 1e-20:
            return None, None
        px = crpix1 + (cd22 * x - cd12 * y) / det
        py = crpix2 + (-cd21 * x + cd11 * y) / det
        return px, py
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return None, None


# ─────────────────────────────────────────────────────────────────────────────
# Catalog queries via VizieR HTTP (no astroquery)
# ─────────────────────────────────────────────────────────────────────────────

def _vizier_tsv(catalog: str, ra: float, dec: float,
                radius_arcmin: float, columns: list[str],
                filters: dict | None = None, max_rows: int = 500) -> list[list[str]]:
    """
    Query VizieR and return parsed data rows (skips comment/header lines).
    Returns [] on any error.
    """
    if not HAS_REQUESTS:
        return []
    params: dict = {
        "-source": catalog,
        "-c": f"{ra:.6f} {dec:+.6f}",
        "-c.r": str(radius_arcmin),
        "-c.u": "arcmin",
        "-out": ",".join(columns),
        "-out.max": str(max_rows),
    }
    if filters:
        params.update(filters)
    try:
        r = requests.get(VIZIER, params=params, timeout=20)
        r.raise_for_status()
        rows = []
        header_skipped = 0
        for line in r.text.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith("-"):
                header_skipped = 0
                continue
            header_skipped += 1
            if header_skipped <= 2:   # column names + units rows
                continue
            rows.append(line.split("\t"))
        return rows
    except Exception:
        return []


def query_vsx(ra: float, dec: float, radius_deg: float) -> list[dict]:
    """Variable stars from VizieR B/vsx/vsx."""
    radius_arcmin = radius_deg * 60
    rows = _vizier_tsv(
        "B/vsx/vsx", ra, dec, radius_arcmin,
        ["Name", "Type", "Period", "max", "RAJ2000", "DEJ2000"],
    )
    stars = []
    for parts in rows:
        if len(parts) < 6:
            continue
        try:
            name     = parts[0].strip()
            var_type = parts[1].strip() or "?"
            period   = parts[2].strip()
            mag_str  = parts[3].strip()
            ra_s     = parts[4].strip()
            dec_s    = parts[5].strip()
            if not ra_s or not dec_s:
                continue
            stars.append({
                "name":     name,
                "ra":       float(ra_s),
                "dec":      float(dec_s),
                "mag":      float(mag_str) if mag_str else 99.0,
                "var_type": var_type,
                "period":   f"{float(period):.4f} d" if period else "—",
                "dist":     0.0,
            })
        except (ValueError, IndexError):
            continue
    return stars


def query_apass(ra: float, dec: float, target_mag: float,
                radius_deg: float = 1.0, n: int = 6) -> list[dict]:
    """Comparison stars from VizieR II/336/apass9."""
    rows = _vizier_tsv(
        "II/336/apass9", ra, dec, radius_deg * 60,
        ["RAJ2000", "DEJ2000", "Vmag", "e_Vmag"],
        filters={"e_Vmag": "<0.05"},
        max_rows=200,
    )
    comps = []
    for parts in rows:
        if len(parts) < 3:
            continue
        try:
            ra_c  = float(parts[0].strip())
            dec_c = float(parts[1].strip())
            vmag  = float(parts[2].strip())
            if abs(vmag - target_mag) > 2.0:
                continue
            # exclude if too close to target (< 5 arcsec)
            sep = ((ra_c - ra) * 3600) ** 2 + ((dec_c - dec) * 3600) ** 2
            if sep < 25:
                continue
            comps.append({"ra": ra_c, "dec": dec_c, "vmag": vmag})
        except (ValueError, IndexError):
            continue
    comps.sort(key=lambda c: abs(c["vmag"] - target_mag))
    return comps[:n]

# ─────────────────────────────────────────────────────────────────────────────
# Siril runner
# ─────────────────────────────────────────────────────────────────────────────

def find_siril_cli() -> str:
    for c in ["siril-cli",
              "/Applications/Siril.app/Contents/MacOS/siril-cli",
              "/usr/local/bin/siril-cli",
              "/opt/homebrew/bin/siril-cli"]:
        if shutil.which(c):
            return c
    return ""


class SirilRunner:
    def __init__(self, log_cb=None):
        self._log = log_cb or print
        self._iface = None
        if HAS_SIRILPY:
            try:
                self._iface = SirilInterface()
                self._iface.connect()
                self._log("Connected to Siril via sirilpy ✓")
            except Exception as e:
                self._log(f"sirilpy: {e} — using siril-cli")

    def run_script(self, working_dir: Path, commands: list[str],
                   name: str = "_phot.ssf") -> bool:
        script = working_dir / name
        script.write_text("requires 1.4\n" + "\n".join(commands) + "\n",
                          encoding="utf-8")
        cli = find_siril_cli()
        if not cli:
            self._log("ERROR: siril-cli not found. Install Siril 1.4+")
            return False
        self._log(f"→ {cli} -d {working_dir.name} -s {name}")
        proc = subprocess.Popen(
            [cli, "-d", str(working_dir), "-s", str(script)],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        )
        for line in proc.stdout:
            self._log(line.rstrip())
        proc.wait()
        return proc.returncode == 0

    def get_working_dir(self) -> Optional[Path]:
        if self._iface:
            try:
                return Path(self._iface.get_wd())
            except Exception:
                pass
        return Path(os.getcwd())


def truncate_comp_csv(csv_path: Path, n: int) -> int:
    """Keep only the first n Comp1 entries in a findcompstars CSV. Returns kept count."""
    try:
        lines = csv_path.read_text(encoding="utf-8").splitlines()
        out, comp_kept = [], 0
        for line in lines:
            if line.startswith("Comp1,") and comp_kept >= n:
                continue
            if line.startswith("Comp1,"):
                comp_kept += 1
            out.append(line)
        csv_path.write_text("\n".join(out) + "\n", encoding="utf-8")
        return comp_kept
    except Exception:
        return 0


def dateobs_to_jd(date_str: str) -> Optional[float]:
    """Convert a FITS DATE-OBS string (ISO 8601) to Julian Date (mid-exposure not applied)."""
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%dT%H:%M"):
        try:
            dt = datetime.datetime.strptime(date_str.strip(), fmt)
            delta = dt - datetime.datetime(2000, 1, 1, 12, 0, 0)
            return 2451545.0 + delta.total_seconds() / 86400.0
        except ValueError:
            continue
    return None


def inject_jd_into_dat(dat_path: Path, seq_path: Path,
                       stem: str, fixlen: int, proc: Path,
                       exptime_s: float = 20.0) -> bool:
    """Replace frame-index JD column in light_curve.dat with real Julian Dates.

    Reads DATE-OBS from each r_light_XXXXX.fit file (selected frames in seq order).
    Adds half the exposure time for mid-exposure JD.
    Returns True if at least one JD was injected.
    """
    # Build ordered list of img_numbers for selected frames
    selected_imgs: list[int] = []
    try:
        with open(seq_path) as fh:
            for line in fh:
                if line.startswith("I "):
                    parts = line.split()
                    img_num, flag = int(parts[1]), int(parts[2])
                    if flag == 1:
                        selected_imgs.append(img_num)
    except Exception:
        return False

    # Map 1-based frame index → JD
    jd_map: dict[int, float] = {}
    half_exp = exptime_s / 86400.0 / 2.0
    for idx, img_num in enumerate(selected_imgs, start=1):
        fit_path = proc / f"{stem}_{img_num:0{fixlen}d}.fit"
        hdr = read_fits_header(fit_path)
        date_str = str(hdr.get("DATE-OBS", ""))
        if date_str:
            jd = dateobs_to_jd(date_str)
            if jd is not None:
                jd_map[idx] = jd + half_exp  # mid-exposure

    if not jd_map:
        return False

    # Rewrite dat file: replace integer frame index with real JD
    try:
        lines = dat_path.read_text(encoding="utf-8").splitlines()
        out: list[str] = []
        for line in lines:
            if "#JD_UT (+ 0)" in line:
                out.append("#JD_UT")
                continue
            if line.startswith("#"):
                out.append(line)
                continue
            parts = line.split()
            if parts:
                try:
                    idx = round(float(parts[0]))
                    if idx in jd_map:
                        parts[0] = f"{jd_map[idx]:.6f}"
                except (ValueError, KeyError):
                    pass
            out.append(" ".join(parts))
        dat_path.write_text("\n".join(out) + "\n", encoding="utf-8")
        return True
    except Exception:
        return False


def get_gain_eadu(lights_dir: Path) -> float:
    """Read real gain in e-/ADU from the first FITS file in lights/.
    GAIN=200 on a Seestar is ISO-equivalent, not e-/ADU.
    EGAIN (if present) is the actual conversion factor.
    Returns 1.0 as a safe fallback for modern CMOS if nothing found."""
    for fit in sorted(lights_dir.glob("*.fit"))[:1]:
        hdr = read_fits_header(fit)
        for key in ("EGAIN", "EPERDN", "GAIN_E", "CCDGAIN"):
            v = hdr.get(key)
            if v is not None:
                try:
                    g = float(v)
                    if 0.05 < g < 30:   # realistic e-/ADU range
                        return round(g, 3)
                except (ValueError, TypeError):
                    pass
        # GAIN keyword: only trust it if it's in e-/ADU range
        v = hdr.get("GAIN")
        if v is not None:
            try:
                g = float(v)
                if 0.05 < g < 30:
                    return round(g, 3)
            except (ValueError, TypeError):
                pass
    return 1.0

# ─────────────────────────────────────────────────────────────────────────────
# Theme
# ─────────────────────────────────────────────────────────────────────────────

BG      = "#16162a"   # main background
BG2     = "#1e1e35"   # widget background
SURFACE = "#2a2a45"   # raised surface
BORDER  = "#3a3a5c"   # borders
FG      = "#e8e8f0"   # primary text  (high contrast on BG)
FG2     = "#b0b0c8"   # secondary text
BLUE    = "#7aa2f7"   # accent / headings
CYAN    = "#7dcfff"   # info values
GREEN   = "#9ece6a"   # success
RED     = "#f7768e"   # error / target star
YELLOW  = "#e0af68"   # warnings / field info
MONO    = ("Menlo", "Courier New", "monospace")

# ─────────────────────────────────────────────────────────────────────────────
# Dark-themed button (tk.Label based)
# macOS Aqua ignores bg/fg on tk.Button entirely — use Label + bindings instead
# ─────────────────────────────────────────────────────────────────────────────

class DarkButton(tk.Label):
    def __init__(self, parent, text, command, hover_bg=None, **kw):
        self._cmd      = command
        self._enabled  = kw.pop("state", "normal") == "normal"
        kw.setdefault("bg",     SURFACE)
        kw.setdefault("fg",     FG)
        kw.setdefault("font",   ("Helvetica", 11))
        kw.setdefault("padx",   10)
        kw.setdefault("pady",   5)
        kw.setdefault("relief", "flat")
        kw.setdefault("cursor", "hand2")
        self._cur_bg   = kw["bg"]
        self._hover_bg = hover_bg or BORDER
        super().__init__(parent, text=text, **kw)
        self.bind("<Button-1>", lambda e: self._enabled and self._cmd())
        self.bind("<Enter>",    lambda e: self._enabled and tk.Label.configure(self, bg=self._hover_bg))
        self.bind("<Leave>",    lambda e: tk.Label.configure(self, bg=self._cur_bg))

    def configure(self, **kw):
        if "state" in kw:
            self._enabled = kw.pop("state") == "normal"
            tk.Label.configure(self, cursor="hand2" if self._enabled else "arrow")
        if "bg" in kw:
            self._cur_bg = kw["bg"]
        tk.Label.configure(self, **kw)

    config = configure

# ─────────────────────────────────────────────────────────────────────────────
# Main window
# ─────────────────────────────────────────────────────────────────────────────

class App(tk.Tk):
    def __init__(self, session_dir: Optional[Path] = None):
        super().__init__()
        self.title("Seestar S30 Pro — Variable Star Photometry")
        self.configure(bg=BG)
        self.minsize(980, 700)

        self.session_dir:   Optional[Path]  = session_dir
        self.field_ra:      float           = 0.0
        self.field_dec:     float           = 0.0
        self.all_stars:     list[dict]      = []
        self.comp_stars:    list[dict]      = []
        self.selected_star: Optional[dict]  = None
        self.runner:        Optional[SirilRunner] = None

        self._build_ui()
        self._apply_ttk_style()
        self.after(100, self._init_bg)

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        # ── Header ──
        hdr = tk.Label(self,
                       text="  Seestar S30 Pro — Variable Star Finder & Photometry",
                       bg=BG, fg=BLUE,
                       font=("Helvetica", 15, "bold"), anchor="w")
        hdr.pack(fill="x", padx=12, pady=(12, 2))

        tk.Frame(self, bg=BORDER, height=1).pack(fill="x", padx=12)

        # ── Session row ──
        row = self._row(self)
        row.pack(fill="x", padx=12, pady=6)
        tk.Label(row, text="Session:", bg=BG, fg=FG2,
                 font=("Helvetica", 12)).pack(side="left")
        self.sess_var = tk.StringVar()
        tk.Entry(row, textvariable=self.sess_var, bg=BG2, fg=FG,
                 insertbackground=FG, font=("Helvetica", 12), width=52,
                 relief="flat", highlightthickness=1,
                 highlightbackground=BORDER,
                 highlightcolor=BLUE).pack(side="left", padx=(6, 4))
        self._button(row, "Browse…",  self._browse).pack(side="left", padx=3)
        self._button(row, "Load",     self._load_session).pack(side="left")

        # ── Field info ──
        info = self._row(self)
        info.pack(fill="x", padx=12, pady=(0, 4))
        self.lbl_field = tk.Label(info, text="Field: —", bg=BG, fg=YELLOW,
                                  font=("Helvetica", 11))
        self.lbl_field.pack(side="left")
        self.lbl_siril = tk.Label(info, text="", bg=BG, fg=FG2,
                                  font=("Helvetica", 10))
        self.lbl_siril.pack(side="right")

        tk.Frame(self, bg=BORDER, height=1).pack(fill="x", padx=12)

        # ── Paned layout ──
        paned = tk.PanedWindow(self, orient="horizontal",
                               bg=BORDER, sashwidth=4, sashrelief="flat")
        paned.pack(fill="both", expand=True, padx=12, pady=8)

        # Left panel — star table
        left = tk.Frame(paned, bg=BG)
        paned.add(left, minsize=500)

        tbl_hdr = self._row(left)
        tbl_hdr.pack(fill="x", pady=(0, 4))
        tk.Label(tbl_hdr, text="Variable stars in field",
                 bg=BG, fg=BLUE, font=("Helvetica", 12, "bold")).pack(side="left")
        self.lbl_count = tk.Label(tbl_hdr, text="0 stars", bg=BG, fg=FG2,
                                  font=("Helvetica", 11))
        self.lbl_count.pack(side="right")
        self._button(tbl_hdr, "↻  Query VSX catalog",
                     self._fetch_vsx).pack(side="right", padx=6)

        tree_frame = tk.Frame(left, bg=BG)
        tree_frame.pack(fill="both", expand=True)
        vsb = ttk.Scrollbar(tree_frame, orient="vertical")
        vsb.pack(side="right", fill="y")
        cols = ("name", "type", "period", "mag", "dist")
        self.tree = ttk.Treeview(tree_frame, columns=cols, show="headings",
                                 selectmode="browse", style="App.Treeview",
                                 yscrollcommand=vsb.set)
        headers = [("name","Name",230), ("type","Type",75),
                   ("period","Period",100), ("mag","Mag",60), ("dist","Dist°",65)]
        for col, label, w in headers:
            self.tree.heading(col, text=label,
                              command=lambda c=col: self._sort(c))
            anchor = "w" if col == "name" else "center"
            self.tree.column(col, width=w, anchor=anchor, stretch=(col=="name"))
        self.tree.pack(side="left", fill="both", expand=True)
        vsb.configure(command=self.tree.yview)
        self.tree.bind("<<TreeviewSelect>>", self._on_select)

        filt_row = self._row(left)
        filt_row.pack(fill="x", pady=(4, 0))
        tk.Label(filt_row, text="Filter:", bg=BG, fg=FG2,
                 font=("Helvetica", 11)).pack(side="left")
        self.filter_var = tk.StringVar()
        self.filter_var.trace_add("write", lambda *_: self._refresh_table())
        tk.Entry(filt_row, textvariable=self.filter_var, bg=BG2, fg=FG,
                 insertbackground=FG, font=("Helvetica", 11), width=18,
                 relief="flat", highlightthickness=1,
                 highlightbackground=BORDER,
                 highlightcolor=BLUE).pack(side="left", padx=6)
        tk.Label(filt_row, text="Max mag:", bg=BG, fg=FG2,
                 font=("Helvetica", 11)).pack(side="left")
        self.mag_max_var = tk.StringVar(value="14")
        self.mag_max_var.trace_add("write", lambda *_: self._refresh_table())
        tk.Entry(filt_row, textvariable=self.mag_max_var, bg=BG2, fg=FG,
                 insertbackground=FG, font=("Helvetica", 11), width=5,
                 relief="flat", highlightthickness=1,
                 highlightbackground=BORDER,
                 highlightcolor=BLUE).pack(side="left", padx=4)

        # Right panel
        right = tk.Frame(paned, bg=BG)
        paned.add(right, minsize=380)

        # Target card
        card = tk.LabelFrame(right, text=" Selected target ",
                             bg=BG, fg=BLUE, font=("Helvetica", 11, "bold"),
                             relief="solid", bd=1, padx=8, pady=6)
        card.pack(fill="x", pady=(0, 8))

        self.lbl_name = tk.Label(card, text="—", bg=BG, fg=RED,
                                 font=("Helvetica", 15, "bold"), anchor="w")
        self.lbl_name.pack(fill="x")
        self.lbl_type = tk.Label(card, text="", bg=BG, fg=CYAN,
                                 font=("Helvetica", 11), anchor="w")
        self.lbl_type.pack(fill="x")
        self.lbl_coord = tk.Label(card, text="", bg=BG, fg=FG2,
                                  font=("Helvetica", 10), anchor="w")
        self.lbl_coord.pack(fill="x")
        tk.Frame(card, bg=BORDER, height=1).pack(fill="x", pady=4)
        self.lbl_comp = tk.Label(card, text="Comparison stars: 0",
                                 bg=BG, fg=FG2, font=("Helvetica", 11), anchor="w")
        self.lbl_comp.pack(fill="x")
        self._button(card, "Fetch comparison stars (APASS)",
                     self._fetch_comp).pack(anchor="w", pady=(4, 0))

        # ── Observation settings ──────────────────────────────────────────────
        obs_card = tk.LabelFrame(right, text=" Observation settings ",
                                 bg=BG, fg=BLUE, font=("Helvetica", 11, "bold"),
                                 relief="solid", bd=1, padx=8, pady=6)
        obs_card.pack(fill="x", pady=(0, 6))

        flt_row = self._row(obs_card)
        flt_row.pack(fill="x", pady=(0, 4))
        tk.Label(flt_row, text="Filter:", bg=BG, fg=FG2,
                 font=("Helvetica", 11), width=9, anchor="w").pack(side="left")
        self.obs_filter_var = tk.StringVar(value="LP anti-pollution filter (CV)")
        ttk.Combobox(flt_row, textvariable=self.obs_filter_var,
                     values=list(FILTER_OPTIONS.keys()),
                     state="readonly", width=30,
                     style="Dark.TCombobox").pack(side="left", padx=(0, 4))

        step_row = self._row(obs_card)
        step_row.pack(fill="x")
        tk.Label(step_row, text="Start from:", bg=BG, fg=FG2,
                 font=("Helvetica", 11), width=9, anchor="w").pack(side="left")
        self.start_from_var = tk.StringVar(value="Full pipeline (steps 1–5)")
        ttk.Combobox(step_row, textvariable=self.start_from_var,
                     values=["Full pipeline (steps 1–5)",
                             "Apply reg + plate solve (steps 3–5)",
                             "Plate solve only (step 4–5)",
                             "Photometry only (step 5)"],
                     state="readonly", width=34,
                     style="Dark.TCombobox").pack(side="left", padx=(0, 4))

        nstars_row = self._row(obs_card)
        nstars_row.pack(fill="x", pady=(4, 0))
        tk.Label(nstars_row, text="Comp stars:", bg=BG, fg=FG2,
                 font=("Helvetica", 11), width=9, anchor="w").pack(side="left")
        self.nstars_var = tk.IntVar(value=10)
        tk.Spinbox(nstars_row, from_=3, to=19, textvariable=self.nstars_var,
                   width=5, bg=BG2, fg=FG, insertbackground=FG,
                   buttonbackground=SURFACE, font=("Helvetica", 11),
                   relief="flat", highlightthickness=1,
                   highlightbackground=BORDER,
                   highlightcolor=BLUE).pack(side="left", padx=(0, 8))
        tk.Label(nstars_row, text="(3–19, Siril max=19)", bg=BG, fg=FG2,
                 font=("Helvetica", 10)).pack(side="left")

        # ── Calibration frames ────────────────────────────────────────────────
        cal_card = tk.LabelFrame(right, text=" Calibration frames  (optional) ",
                                 bg=BG, fg=BLUE, font=("Helvetica", 11, "bold"),
                                 relief="solid", bd=1, padx=8, pady=6)
        cal_card.pack(fill="x", pady=(0, 8))

        self.dark_var = tk.StringVar()
        self.flat_var = tk.StringVar()
        self.bias_var = tk.StringVar()
        for label, var in [("Darks:", self.dark_var),
                            ("Flats:", self.flat_var),
                            ("Bias:",  self.bias_var)]:
            row = self._row(cal_card)
            row.pack(fill="x", pady=2)
            tk.Label(row, text=label, bg=BG, fg=FG2,
                     font=("Helvetica", 11), width=7, anchor="w").pack(side="left")
            tk.Entry(row, textvariable=var, bg=BG2, fg=FG,
                     insertbackground=FG, font=("Helvetica", 10), width=22,
                     relief="flat", highlightthickness=1,
                     highlightbackground=BORDER,
                     highlightcolor=BLUE).pack(side="left", padx=(0, 4))
            self._button(row, "Browse…",
                         lambda v=var: self._browse_calib(v)).pack(side="left")

        # Run button — starts disabled (dark), turns bright green when a target is selected
        self.run_btn = DarkButton(
            right,
            text="▶   Prepare frames & Generate light curve",
            command=self._run_pipeline,
            state="disabled",
            bg=SURFACE,
            fg=FG2,
            hover_bg="#16a34a",
            font=("Helvetica", 13, "bold"),
            pady=10,
        )
        self.run_btn.pack(fill="x", pady=8)

        # Progress
        self.prog_var = tk.IntVar()
        prog = ttk.Progressbar(right, variable=self.prog_var, maximum=100,
                               style="Green.Horizontal.TProgressbar")
        prog.pack(fill="x", pady=(0, 3))
        self.lbl_prog = tk.Label(right, text="Ready", bg=BG, fg=GREEN,
                                 font=("Helvetica", 11), anchor="w")
        self.lbl_prog.pack(fill="x")

        # Siril log — normal state so text is selectable and copyable
        log_card = tk.LabelFrame(right, text=" Siril log  (select & copy with Cmd+C) ",
                                 bg=BG, fg=BLUE, font=("Helvetica", 11, "bold"),
                                 relief="solid", bd=1)
        log_card.pack(fill="both", expand=True, pady=(8, 0))

        log_frame = tk.Frame(log_card, bg=BG2)
        log_frame.pack(fill="both", expand=True, padx=4, pady=4)
        log_sb = ttk.Scrollbar(log_frame, orient="vertical")
        log_sb.pack(side="right", fill="y")
        self.log = tk.Text(log_frame, bg=BG2, fg=FG2,
                           font=(MONO[0], 10),
                           relief="flat", wrap="none",
                           state="normal",     # stays normal → fully selectable
                           height=10, yscrollcommand=log_sb.set)
        self.log.pack(side="left", fill="both", expand=True)
        log_sb.configure(command=self.log.yview)

    def _apply_ttk_style(self):
        s = ttk.Style()
        s.theme_use("clam")
        s.configure("App.Treeview",
                    background=BG2, foreground=FG,
                    fieldbackground=BG2, rowheight=23,
                    font=("Helvetica", 11))
        s.configure("App.Treeview.Heading",
                    background=SURFACE, foreground=BLUE,
                    font=("Helvetica", 11, "bold"), relief="flat")
        s.map("App.Treeview",
              background=[("selected", SURFACE)],
              foreground=[("selected", BLUE)])
        s.configure("Green.Horizontal.TProgressbar",
                    troughcolor=BG2, background="#22c55e",
                    lightcolor="#16a34a", darkcolor="#15803d")
        s.configure("TScrollbar",
                    background=SURFACE, troughcolor=BG2,
                    arrowcolor=FG2, bordercolor=BORDER)
        s.configure("Dark.TCombobox",
                    fieldbackground=BG2, background=SURFACE,
                    foreground=FG, arrowcolor=FG2,
                    bordercolor=BORDER, lightcolor=BG, darkcolor=BG)
        s.map("Dark.TCombobox",
              fieldbackground=[("readonly", BG2)],
              selectbackground=[("readonly", SURFACE)],
              selectforeground=[("readonly", FG)])

    def _row(self, parent):
        return tk.Frame(parent, bg=BG)

    def _button(self, parent, text, cmd):
        return DarkButton(parent, text=text, command=cmd)

    def _browse_calib(self, var: tk.StringVar):
        d = filedialog.askdirectory(title="Select folder containing calibration frames")
        if d:
            var.set(d)

    # ── Init ──────────────────────────────────────────────────────────────────

    def _init_bg(self):
        threading.Thread(target=self._connect_and_detect, daemon=True).start()

    def _connect_and_detect(self):
        self.runner = SirilRunner(log_cb=self._log)
        cli = find_siril_cli()
        if HAS_SIRILPY and self.runner._iface:
            status, color = "● Siril connected", GREEN
        elif cli:
            status, color = "● siril-cli ready", CYAN
        else:
            status, color = "✗ siril-cli not found", RED
        self.after(0, lambda: self.lbl_siril.configure(text=status, fg=color))

        # Auto-detect session
        if self.session_dir is None:
            wd = self.runner.get_working_dir()
            if wd and (wd / "lights").is_dir():
                self.session_dir = wd
            else:
                for base in [Path("/Volumes/Seestar/MyWorks"),
                             Path.home() / "Desktop",
                             Path.home() / "Documents"]:
                    if not base.exists():
                        continue
                    for child in sorted(base.iterdir(), reverse=True):
                        if child.is_dir() and (child / "lights").is_dir():
                            self.session_dir = child
                            break
                    if self.session_dir:
                        break

        if self.session_dir:
            self.after(0, lambda: self.sess_var.set(str(self.session_dir)))
            self.after(0, self._load_session)

    # ── Session ───────────────────────────────────────────────────────────────

    def _browse(self):
        d = filedialog.askdirectory(title="Select session folder",
                                    initialdir=str(Path.home()))
        if d:
            self.sess_var.set(d)
            self._load_session()

    def _load_session(self):
        s = self.sess_var.get().strip()
        if not s:
            return
        self.session_dir = Path(s)
        if not self.session_dir.is_dir():
            messagebox.showerror("Not found", str(self.session_dir))
            return
        threading.Thread(target=self._load_session_bg, daemon=True).start()

    def _load_session_bg(self):
        session = self.session_dir
        stack = None
        for cand in ["process/lights.fit", "process/result.fit"]:
            p = session / cand
            if p.exists():
                stack = p
                break
        if stack is None:
            for p in session.glob("*_og.fit"):
                stack = p
                break

        ra = dec = 0.0
        obj = ""
        n_frames = 0
        if stack:
            h = read_fits_header(stack)
            ra       = float(h.get("RA",       h.get("CRVAL1", 0)))
            dec      = float(h.get("DEC",      h.get("CRVAL2", 0)))
            obj      = str(h.get("OBJECT",     "?"))
            n_frames = int(float(h.get("STACKCNT", 0)))

        self.field_ra, self.field_dec = ra, dec
        info = f"{obj}  ·  RA {ra:.4f}°  Dec {dec:+.4f}°  ·  {n_frames} stacked frames"
        self.after(0, lambda: self.lbl_field.configure(text=f"Field: {info}"))

        # Load starsv.csv
        stars: list[dict] = []
        csv_path = session / "starsv.csv"
        if csv_path.exists():
            with open(csv_path, newline="", encoding="utf-8") as f:
                for row in csv.DictReader(f):
                    try:
                        sr = float(row.get("ra", 0))
                        sd = float(row.get("dec", 0))
                        dist = 0.0
                        if ra:
                            dist = round(
                                (((sr - ra) * 0.7071) ** 2 + (sd - dec) ** 2) ** 0.5, 3
                            )
                        stars.append({
                            "name":     row.get("name", "?"),
                            "ra":       sr, "dec": sd,
                            "mag":      float(row.get("mag", 99) or 99),
                            "var_type": "?", "period": "—", "dist": dist,
                        })
                    except (ValueError, KeyError):
                        continue
            self._log(f"{len(stars)} stars loaded from starsv.csv")
        else:
            self._log("starsv.csv not found — click '↻ Query VSX catalog'")

        self.all_stars = stars
        self.after(0, self._refresh_table)

        # Auto-select best start point based on what's already done
        proc = session / "process"
        has_registered = any(
            (proc / f"{seq}.seq").exists()
            for seq in ("r_light_", "r_pp_light_")
        )
        has_platesolved = False
        if has_registered:
            for seq in ("r_light_", "r_pp_light_"):
                fits_list = sorted(proc.glob(f"{seq.rstrip('_')}*.fit"))
                if fits_list:
                    try:
                        hdr = read_fits_header(fits_list[0])
                        has_platesolved = "CRVAL1" in hdr
                    except Exception:
                        pass
                    break

        has_lc = (proc / "light_curve.dat").exists()
        if has_platesolved and has_lc:
            suggested = "Photometry only (step 5)"
        elif has_platesolved and not has_lc:
            # Plate solve done but photometry failed — most likely a framing/drift
            # issue in seqapplyreg; re-run from step 3 with -framing=cog
            suggested = "Apply reg + plate solve (steps 3–5)"
        elif has_registered:
            suggested = "Plate solve only (step 4–5)"
        else:
            suggested = "Full pipeline (steps 1–5)"

        self.after(0, lambda s=suggested: self.start_from_var.set(s))

    # ── Table ─────────────────────────────────────────────────────────────────

    def _refresh_table(self):
        f = self.filter_var.get().lower()
        try:
            mag_max = float(self.mag_max_var.get())
        except (ValueError, AttributeError):
            mag_max = 99.0
        self.tree.delete(*self.tree.get_children())
        shown = 0
        for s in self.all_stars:
            if f and f not in s["name"].lower() and f not in s["var_type"].lower():
                continue
            if s.get("mag", 99.0) > mag_max:
                continue
            tag = "known" if s.get("var_type", "?") not in ("?", "", " ") else ""
            self.tree.insert("", "end", iid=str(id(s)), tags=(tag,), values=(
                s["name"],
                s.get("var_type", "?"),
                s.get("period", "—"),
                f"{s['mag']:.2f}",
                f"{s.get('dist', 0):.2f}",
            ))
            shown += 1
        self.tree.tag_configure("known", foreground=CYAN)
        self.lbl_count.configure(text=f"{shown} stars")

    def _sort(self, col: str):
        key_map = {"name": "name", "type": "var_type",
                   "period": "period", "mag": "mag", "dist": "dist"}
        key = key_map.get(col, col)
        rev = getattr(self, f"_rev_{col}", False)
        try:
            self.all_stars.sort(
                key=lambda s: float(str(s.get(key, 0)).split()[0] or 0), reverse=rev
            )
        except (ValueError, TypeError):
            self.all_stars.sort(key=lambda s: str(s.get(key, "")), reverse=rev)
        setattr(self, f"_rev_{col}", not rev)
        self._refresh_table()

    def _on_select(self, _=None):
        sel = self.tree.selection()
        if not sel:
            return
        star = next((s for s in self.all_stars if str(id(s)) == sel[0]), None)
        if not star:
            return
        self.selected_star = star
        self.comp_stars = []
        self.lbl_name.configure(text=star["name"])
        self.lbl_type.configure(
            text=f"Type: {star.get('var_type','?')}   Period: {star.get('period','—')}"
        )
        self.lbl_coord.configure(
            text=f"RA {star['ra']:.5f}°  /  Dec {star['dec']:+.5f}°  ·  mag {star['mag']:.2f}"
        )
        self.lbl_comp.configure(text="Comparison stars: 0")
        self.run_btn.configure(state="normal", bg="#22c55e", fg="#0a2a0a")

    # ── VSX catalog ───────────────────────────────────────────────────────────

    def _fetch_vsx(self):
        if not self.field_ra:
            messagebox.showinfo("Info", "Load a session first.")
            return
        if not HAS_REQUESTS:
            messagebox.showerror("Missing",
                                 "requests not available.\n"
                                 "Run: test_connections.py for diagnostics.")
            return
        self._log("Querying VizieR VSX catalog…")
        self._btn_state(False)
        threading.Thread(target=self._fetch_vsx_bg, daemon=True).start()

    def _fetch_vsx_bg(self):
        radius = max(SEESTAR["fov_w"], SEESTAR["fov_h"]) * 0.75
        stars = query_vsx(self.field_ra, self.field_dec, radius)
        if not stars:
            self.after(0, lambda: (
                self._log("VSX: no results — check connection with test_connections.py"),
                self._btn_state(True),
            ))
            return
        for s in stars:
            s["dist"] = round(
                (((s["ra"] - self.field_ra) * 0.7071) ** 2
                 + (s["dec"] - self.field_dec) ** 2) ** 0.5, 3
            )
        self.all_stars = sorted(stars, key=lambda s: s["dist"])
        self.after(0, lambda: (
            self._refresh_table(),
            self._log(f"VSX: {len(stars)} variable stars found in field"),
            self._btn_state(True),
        ))

    # ── Comparison stars ──────────────────────────────────────────────────────

    def _fetch_comp(self):
        if not self.selected_star:
            messagebox.showinfo("Info", "Select a target star first.")
            return
        self._log(f"Fetching comparison stars for {self.selected_star['name']} (APASS)…")
        threading.Thread(target=self._fetch_comp_bg, daemon=True).start()

    def _fetch_comp_bg(self):
        s = self.selected_star
        comps = query_apass(s["ra"], s["dec"], s["mag"])
        self.comp_stars = comps
        msg = f"{len(comps)} comparison stars from APASS"
        self.after(0, lambda: (
            self.lbl_comp.configure(text=msg),
            self._log(msg),
        ))

    # ── Pipeline ──────────────────────────────────────────────────────────────

    def _run_pipeline(self):
        if not self.selected_star:
            return
        self.run_btn.configure(state="disabled", bg=SURFACE, fg=FG2)
        self.prog_var.set(0)
        self.log.delete("1.0", "end")
        threading.Thread(target=self._pipeline_bg, daemon=True).start()

    def _pipeline_bg(self):
        session = self.session_dir
        proc    = session / "process"
        lights  = session / "lights"
        masters = proc / "masters"
        proc.mkdir(exist_ok=True)

        star  = self.selected_star
        focal = SEESTAR["focal"]
        pixsz = SEESTAR["pixsz"]
        gain  = get_gain_eadu(lights)
        self._log(f"Gain: {gain} e-/ADU")

        # Determine which steps to run
        start_label = self.start_from_var.get()
        if "step 5" in start_label or "Photometry" in start_label:
            start_step = 5
        elif "step 4" in start_label or "Plate solve" in start_label:
            start_step = 4
        elif "step 3" in start_label or "Apply reg" in start_label:
            start_step = 3
        else:
            start_step = 1

        # ─── Jump to step 3/4/5 — detect existing registered sequence ──────
        if start_step == 3:
            # Need light_.seq from the register step
            for candidate in ("pp_light_", "light_"):
                if (proc / f"{candidate}.seq").exists():
                    seq = candidate
                    break
            else:
                self._done(False,
                    f"No registered sequence (light_.seq) found in {proc.name}/.\n"
                    "Run the full pipeline first (steps 1–5)."); return
            registered = f"r_{seq}"

        if start_step >= 4:
            for candidate in ("r_pp_light_", "r_light_"):
                if (proc / f"{candidate}.seq").exists():
                    registered = candidate
                    break
            else:
                self._done(False,
                    f"No registered sequence found in {proc.name}/.\n"
                    "Run the full pipeline first (steps 1–5)."); return
            self._log(f"Resuming from step {start_step} — sequence: {registered}")

        # ─── Steps 0–3 (only when starting from step 1) ──────────────────
        if start_step == 1:
            dark_dir = Path(self.dark_var.get()) if self.dark_var.get().strip() else None
            flat_dir = Path(self.flat_var.get()) if self.flat_var.get().strip() else None
            bias_dir = Path(self.bias_var.get()) if self.bias_var.get().strip() else None
            has_calib = any([dark_dir, flat_dir, bias_dir])

            # Cleanup: remove misplaced directory from previous run with quoted -out
            bad_out = lights / f'"{proc}"'
            if bad_out.exists():
                shutil.rmtree(bad_out, ignore_errors=True)
                self._log("Cleaned up misplaced output directory from a previous run.")

            # Step 0 (optional): stack calibration masters + calibrate lights
            if has_calib:
                masters.mkdir(exist_ok=True)
                self._prog(5, "Building calibration masters…")
                self._log("─── Step 0: calibration frames ───")
                calib_cmds = []
                for kind, src in [("dark", dark_dir), ("flat", flat_dir), ("bias", bias_dir)]:
                    if src is None:
                        continue
                    norm = "-nonorm" if kind in ("dark", "bias") else "-norm=mul"
                    calib_cmds += [
                        f'cd "{src}"',
                        f'link {kind} -out={masters}',  # no quotes: Siril takes them literally
                        f'cd "{masters}"',
                        f'stack {kind}_ rej 3 3 {norm} -out=master_{kind}',
                    ]
                    self._log(f"  → stacking {kind}s from {src.name}/")
                cal_flags = " ".join(
                    f"-{k}=masters/master_{k}"
                    for k, d in [("bias", bias_dir), ("dark", dark_dir), ("flat", flat_dir)]
                    if d is not None
                )
                calib_cmds += [f'cd "{proc}"', f"calibrate light_ {cal_flags} -cc=banding"]
                ok = self.runner.run_script(proc, calib_cmds, "_s0_calibrate.ssf")
                if not ok:
                    self._done(False, "Step 0 failed (calibration)"); return
                seq = "pp_light_"
            else:
                self._log("No calibration frames — proceeding with raw lights.")
                seq = "light_"

            # Steps 1+2: link + register in the same siril-cli session
            # IMPORTANT: -out= path must NOT be quoted — Siril parses quotes literally
            self._prog(15, "Linking frames & computing registration…")
            self._log("─── Step 1: link lights → sequence ───")
            self._log("─── Step 2: register -2pass ───")
            ok = self.runner.run_script(proc, [
                f'cd "{lights}"',
                f'link light -out={proc}',
                f'cd "{proc}"',
                f"register {seq} -2pass",
            ], "_s12_link_register.ssf")
            if not ok:
                self._done(False, "Step 1/2 failed (link / register)"); return

            # Step 3: apply registration
            self._prog(45, "Aligning frames…")
            self._log("─── Step 3: seqapplyreg ───")
            ok = self.runner.run_script(proc, [
                f'cd "{proc}"',
                # -framing=cog: centers output on gravity center of all frames,
                # keeping stored registration shifts small → avoids light_curve
                # "heavy drift" failure that occurs with -framing=max
                f"seqapplyreg {seq} -framing=cog -filter-round=2.5k",
            ], "_s3_applyreg.ssf")
            if not ok:
                self._done(False, "Step 3 failed (seqapplyreg)"); return

            registered = f"r_{seq}"

        if start_step == 3:
            # Came here from "Apply reg + plate solve (steps 3–5)"
            self._prog(45, "Aligning frames…")
            self._log("─── Step 3: seqapplyreg ───")
            ok = self.runner.run_script(proc, [
                f'cd "{proc}"',
                f"seqapplyreg {seq} -framing=cog -filter-round=2.5k",
            ], "_s3_applyreg.ssf")
            if not ok:
                self._done(False, "Step 3 failed (seqapplyreg)"); return

        # ─── Step 4: plate solve (skippable) ─────────────────────────────
        if start_step <= 4:
            self._prog(65, "Plate solving registered frames…")
            self._log("─── Step 4: seqplatesolve ───")
            disto = proc / "ps_distortion"
            disto_arg = "-disto=ps_distortion" if disto.is_dir() else ""
            ok = self.runner.run_script(proc, [
                f'cd "{proc}"',
                f"seqplatesolve {registered} -nocache -force "
                f"-focal={focal} -pixelsize={pixsz} -radius=2.5 {disto_arg}".strip(),
            ], "_s4_platesolve.ssf")
            if not ok:
                self._log("WARNING: plate solve had errors — continuing")

        # ── Step 5: clean stale files, resolve reference frame ───────────────
        self._prog(82, "Checking sequence integrity…")

        # Parse the .seq S-line to find the reference frame image number.
        # S-line: S 'name' start nb_images nb_selected fixed_len ref_idx version ...
        # ref_idx is a 0-based index into the I-lines list.
        seq_file  = proc / f"{registered}.seq"
        stem      = registered.rstrip("_")
        seq_fixlen  = 4
        seq_ref_idx = None          # 0-based index in I-lines
        img_entries: list[tuple[int, int]] = []   # (img_num, 1-based rank)
        try:
            with open(seq_file) as fh:
                rank = 0
                for line in fh:
                    if line.startswith("S "):
                        parts = line.split()
                        seq_fixlen  = int(parts[5])
                        seq_ref_idx = int(parts[6])
                    elif line.startswith("I "):
                        rank += 1
                        img_entries.append((int(line.split()[1]), rank))
        except Exception:
            pass

        # Remove stale .fit files no longer in the .seq
        if img_entries:
            valid_nums = {n for n, _ in img_entries}
            for f in sorted(proc.glob(f"{stem}_*.fit")):
                try:
                    num = int(f.stem[len(stem)+1:])
                    if num not in valid_nums:
                        self._log(f"Removing stale frame {f.name}")
                        f.unlink()
                except (ValueError, OSError):
                    pass

        ref_img_num = None
        if seq_ref_idx is not None and seq_ref_idx < len(img_entries):
            ref_img_num = img_entries[seq_ref_idx][0]

        self._prog(85, "Aperture photometry…")
        self._log("─── Step 5: setphot + light_curve ───")

        comp_csv = proc / "comp_stars.csv"
        lc_cmd   = None

        # Primary: findcompstars → -ninastars -autoring.
        # findcompstars queries APASS catalog, filters variable stars from GCVS,
        # and writes a NINA-format CSV. light_curve -ninastars handles all
        # coordinate conversion internally via the sequence's WCS.
        # Requires a plate-solved reference frame loaded into gfit.
        ref_fits_name = (
            f"{stem}_{ref_img_num:0{seq_fixlen}d}.fit"
            if ref_img_num is not None else None
        )
        if ref_fits_name and (proc / ref_fits_name).exists():
            if comp_csv.exists():
                comp_csv.unlink()
            star_arg = star["name"].replace('"', '\\"')
            self._log(f"findcompstars: querying APASS for {star['name']}…")
            self.runner.run_script(proc, [
                f'cd "{proc}"',
                f"load {ref_fits_name}",
                f'findcompstars "{star_arg}" -narrow -dvmag=3 -emag=0.05 -catalog=apass -out=comp_stars.csv',
            ], "_s5a_findcomp.ssf")
            if comp_csv.exists() and comp_csv.stat().st_size > 50:
                n = max(3, min(50, self.nstars_var.get()))
                kept = truncate_comp_csv(comp_csv, n)
                lc_cmd = f"light_curve {registered} 0 -ninastars=comp_stars.csv"
                self._log(f"Siril comparison stars ready: {kept} stars (findcompstars)")
            else:
                self._log("findcompstars produced no output — falling back to manual comp stars")

        # Fallback A: manual display-space pixel coords computed from our VizieR APASS query.
        # -at/-refat expect Siril display coords: display_x = fits_x - 0.5,
        # display_y = NAXIS2 - fits_y + 0.5  (Y-flip + 0.5 offset, integer-rounded).
        # -autoring is incompatible with -at mode in Siril 1.4.3.
        if lc_cmd is None and ref_img_num is not None and self.comp_stars:
            ref_fits = proc / f"{stem}_{ref_img_num:0{seq_fixlen}d}.fit"
            ref_hdr  = read_fits_header(ref_fits)
            naxis2   = int(ref_hdr.get("NAXIS2", 0))
            tx, ty   = sky_to_pixel(star["ra"], star["dec"], ref_hdr)
            if tx is not None and naxis2 > 0:
                def _to_disp(fx, fy, n2=naxis2):
                    return round(fx - 0.5), round(n2 - fy + 0.5)
                ref_pix = []
                for cs in self.comp_stars:
                    rx, ry = sky_to_pixel(cs["ra"], cs["dec"], ref_hdr)
                    if rx is not None:
                        ref_pix.append(_to_disp(rx, ry))
                if ref_pix:
                    tdx, tdy = _to_disp(tx, ty)
                    lc_cmd = f"light_curve {registered} 0 -at={tdx},{tdy}"
                    for rdx, rdy in ref_pix:
                        lc_cmd += f" -refat={rdx},{rdy}"
                    self._log(
                        f"Fallback: display coords target ({tdx},{tdy}), "
                        f"{len(ref_pix)} comp stars"
                    )

        # Fallback B: sky coordinates (-wcs/-refwcs).
        if lc_cmd is None:
            self._log("WCS failed — using sky coordinates (-wcs)")
            lc_cmd = (f"light_curve {registered} 0 "
                      f"-wcs={star['ra']:.6f},{star['dec']:.6f}")
            for cs in self.comp_stars:
                lc_cmd += f" -refwcs={cs['ra']:.6f},{cs['dec']:.6f}"

        ok = self.runner.run_script(proc, [
            f'cd "{proc}"',
            f"setphot -aperture=8 -inner=14 -outer=21 -dyn_ratio=4.0 -gain={gain}",
            lc_cmd,
        ], "_s5_phot.ssf")
        if not ok:
            self._done(False, "Step 5 failed (light_curve)"); return

        # ── Collect results ──
        lc_dat = proc / "light_curve.dat"
        results_dir = session / "results"
        results_dir.mkdir(exist_ok=True)
        safe = star["name"].replace(" ", "_").replace("/", "-")
        out_dat = results_dir / f"{safe}_light_curve.dat"
        out_csv = results_dir / f"{safe}_aavso.csv"

        if lc_dat.exists():
            # Inject real Julian Dates from DATE-OBS headers before copying
            if inject_jd_into_dat(lc_dat, seq_file, stem, seq_fixlen, proc):
                self._log("JD timestamps injected from DATE-OBS headers ✓")
            else:
                self._log("WARNING: DATE-OBS not found — JD axis shows frame indices")
            shutil.copy2(lc_dat, out_dat)
            png = proc / "light_curve.png"
            if png.exists():
                shutil.copy2(png, results_dir / f"{safe}_light_curve.png")
            self._export_aavso(out_dat, out_csv, star["name"])
            self._done(True, str(out_dat))
        else:
            self._done(False,
                       "light_curve.dat not found.\n"
                       "The target may be outside the field, or too faint.")

    def _export_aavso(self, dat: Path, out: Path, name: str):
        rows = []
        with open(dat) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                try:
                    vals = [float(p) for p in parts]
                    jd_i = next(i for i, v in enumerate(vals) if v > 2_400_000)
                    rows.append((vals[jd_i], vals[jd_i+1],
                                 vals[jd_i+2] if len(vals) > jd_i+2 else 0.0))
                except (ValueError, StopIteration):
                    continue
        if not rows:
            return
        filt_key  = self.obs_filter_var.get()
        filt_code, filt_note = FILTER_OPTIONS.get(filt_key, ("CV", ""))
        notes = f"seestar_s30pro|{filt_note}" if filt_note else "seestar_s30pro"

        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            for d in ["#TYPE=EXTENDED", "#OBSCODE=XXXX",
                      "#SOFTWARE=Siril+seestar_varstar_siril.py",
                      f"#FILTER={filt_code}",
                      "#DELIM=,", "#DATE=JD", "#OBSTYPE=CCD"]:
                w.writerow([d])
            w.writerow(["NAME","DATE","MAG","MERR","FILT","TRANS","MTYPE",
                        "CNAME","CMAG","KNAME","KMAG","AMASS","GROUP","CHART","NOTES"])
            for jd, mag, err in rows:
                # MTYPE=DIFF: differential photometry — no instrumental→standard transform
                w.writerow([name, f"{jd:.6f}", f"{mag:.4f}", f"{err:.4f}",
                            filt_code, "NO", "DIFF",
                            "ENSEMBLE", "na", "na", "na", "na", "1", "na",
                            notes])

    # ── Thread-safe helpers ───────────────────────────────────────────────────

    def _log(self, msg: str):
        def _do():
            self.log.insert("end", msg + "\n")
            self.log.see("end")
        self.after(0, _do)

    def _prog(self, pct: int, label: str):
        self.after(0, lambda: (
            self.prog_var.set(pct),
            self.lbl_prog.configure(text=label, fg=CYAN),
        ))

    def _btn_state(self, enabled: bool):
        if self.selected_star:
            self.after(0, lambda: self.run_btn.configure(
                state="normal" if enabled else "disabled",
                bg="#22c55e" if enabled else SURFACE,
                fg="#0a2a0a" if enabled else FG2,
            ))

    def _done(self, success: bool, msg: str):
        self.after(0, lambda: self._done_ui(success, msg))

    def _done_ui(self, ok: bool, msg: str):
        self.run_btn.configure(state="normal", bg="#22c55e", fg="#0a2a0a")
        self.prog_var.set(100 if ok else 0)
        if ok:
            self.lbl_prog.configure(text="Done!", fg=GREEN)
            self._log(f"✓  {msg}")
            messagebox.showinfo("Done!",
                                f"Light curve saved:\n{msg}\n\n"
                                "Results are in the results/ folder.")
        else:
            self.lbl_prog.configure(text=f"Failed: {msg}", fg=RED)
            self._log(f"✗  {msg}")
            messagebox.showerror("Pipeline error", msg)

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--session", help="Session folder (auto-detected if omitted)")
    args, _ = p.parse_known_args()
    session = Path(args.session) if args.session else None
    App(session_dir=session).mainloop()
