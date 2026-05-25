#!/usr/bin/env python3
"""
Seestar S30 Pro — Variable Star Finder & Photometry GUI
========================================================
A PyQt6 application to find variable stars in a Seestar session field and
generate differential light curves for a selected star using Siril.

Run it from terminal (Siril must be installed with siril-cli in PATH):
    python3 seestar_varstar_gui.py
    python3 seestar_varstar_gui.py --session /Volumes/Seestar/MyWorks/postprod_m81

Requirements:
    pip3 install PyQt6 astropy astroquery
    Siril 1.4+ with siril-cli in PATH
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

# ── PyQt6 ────────────────────────────────────────────────────────────────────
try:
    from PyQt6.QtCore import (
        Qt, QThread, pyqtSignal, QAbstractTableModel, QModelIndex, QSortFilterProxyModel,
    )
    from PyQt6.QtGui import QColor, QFont, QPalette, QIcon
    from PyQt6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
        QLabel, QLineEdit, QPushButton, QTableView, QHeaderView, QTextEdit,
        QProgressBar, QGroupBox, QComboBox, QFileDialog, QMessageBox, QSplitter,
        QStatusBar, QFrame, QSpinBox, QCheckBox, QTabWidget,
    )
except ImportError:
    print("ERROR: PyQt6 not installed. Run: pip3 install PyQt6")
    sys.exit(1)

# ── Science deps ──────────────────────────────────────────────────────────────
try:
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False
    print("WARNING: astropy not installed. pip3 install astropy")

try:
    from astroquery.vizier import Vizier
    from astroquery.simbad import Simbad
    HAS_ASTROQUERY = True
except ImportError:
    HAS_ASTROQUERY = False
    print("WARNING: astroquery not installed. pip3 install astroquery")

try:
    import matplotlib
    matplotlib.use("QtAgg")
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
    from matplotlib.figure import Figure
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# ─────────────────────────────────────────────────────────────────────────────
# Hardware profile — updated from actual FITS headers
# ─────────────────────────────────────────────────────────────────────────────
SEESTAR_PROFILE = {
    "focal_mm": 160,       # from FITS FOCALLEN header
    "pixel_um": 2.9,       # from FITS XPIXSZ header
    "gain": 200,           # from FITS GAIN header
    "sensor": "IMX585",
    "bayer": "GRBG",
    "filter": "LP",
    "fov_deg_w": 2.2,      # ~2160px × 2.9µm / 160mm in degrees
    "fov_deg_h": 3.9,      # ~3840px × 2.9µm / 160mm in degrees
}

# Default photometry channel: green (best SNR on GRBG Bayer OSC)
PHOT_CHANNEL = 1
PHOT_APERTURE = 10
PHOT_INNER = 15
PHOT_OUTER = 25

# Auto-scan paths for Seestar USB
SEESTAR_AUTO_PATHS = [
    "/Volumes/Seestar/MyWorks",
    os.path.expanduser("~/Desktop"),
    os.path.expanduser("~/Documents"),
]


# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class FieldInfo:
    ra: float = 0.0
    dec: float = 0.0
    object_name: str = ""
    n_frames: int = 0
    exptime: float = 0.0
    date_obs: str = ""
    stack_path: Optional[Path] = None
    wcs: Optional[object] = None

    @property
    def coord(self) -> Optional[object]:
        if not HAS_ASTROPY:
            return None
        return SkyCoord(self.ra, self.dec, unit=u.deg)

    @property
    def ra_hms(self) -> str:
        if not HAS_ASTROPY:
            return f"{self.ra:.4f}°"
        return self.coord.ra.to_string(unit=u.hour, sep=":", precision=1)

    @property
    def dec_dms(self) -> str:
        if not HAS_ASTROPY:
            return f"{self.dec:+.4f}°"
        return self.coord.dec.to_string(sep=":", precision=0, alwayssign=True)


@dataclass
class VariableStar:
    name: str
    ra: float
    dec: float
    mag: float
    var_type: str = "?"
    period: str = "—"
    amplitude: str = "—"
    in_fov: bool = True
    dist_from_center: float = 0.0  # degrees

    def wcs_arg(self) -> str:
        return f"{self.ra:.6f},{self.dec:+.6f}"


@dataclass
class CompStar:
    ra: float
    dec: float
    vmag: float
    name: str = ""


# ─────────────────────────────────────────────────────────────────────────────
# Siril CLI runner
# ─────────────────────────────────────────────────────────────────────────────

def find_siril_cli() -> str:
    candidates = [
        "siril-cli",
        "/usr/local/bin/siril-cli",
        "/opt/homebrew/bin/siril-cli",
        "/Applications/Siril.app/Contents/MacOS/siril-cli",
        "/usr/bin/siril-cli",
    ]
    for c in candidates:
        if shutil.which(c):
            return c
    return ""


def run_siril(working_dir: Path, commands: List[str],
              log_callback=None, script_name: str = "_varphot.ssf") -> bool:
    """Execute a Siril .ssf script. Streams output lines to log_callback."""
    script_path = working_dir / script_name
    script_path.write_text("requires 1.4\n" + "\n".join(commands) + "\n", encoding="utf-8")

    cli = find_siril_cli()
    if not cli:
        if log_callback:
            log_callback("ERROR: siril-cli not found. Install Siril 1.4+.")
        return False

    cmd = [cli, "-d", str(working_dir), "-s", str(script_path)]
    if log_callback:
        log_callback(f"→ {' '.join(cmd)}")

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    for line in proc.stdout:
        line = line.rstrip()
        if line and log_callback:
            log_callback(line)
    proc.wait()
    return proc.returncode == 0


# ─────────────────────────────────────────────────────────────────────────────
# Worker threads
# ─────────────────────────────────────────────────────────────────────────────

class SessionLoader(QThread):
    """Reads stacked FITS and existing starsv.csv from a session folder."""
    done = pyqtSignal(object, list)   # FieldInfo, [VariableStar]
    error = pyqtSignal(str)

    def __init__(self, session_dir: Path):
        super().__init__()
        self.session_dir = session_dir

    def run(self):
        try:
            fi, stars = _load_session(self.session_dir)
            self.done.emit(fi, stars)
        except Exception as exc:
            self.error.emit(str(exc))


def _load_session(session_dir: Path) -> Tuple[FieldInfo, List[VariableStar]]:
    """Parse session folder: read WCS from stack, load starsv.csv."""
    fi = FieldInfo()

    # Find stacked image (lights.fit, result.fit, *_og.fit …)
    stack = None
    for candidate in ["process/lights.fit", "process/result.fit"]:
        p = session_dir / candidate
        if p.exists():
            stack = p
            break
    if stack is None:
        # Try any .fit in process/
        fits_in_process = sorted((session_dir / "process").glob("*.fit"))
        for p in fits_in_process:
            if "light" in p.stem.lower() or "result" in p.stem.lower():
                stack = p
                break
        if not stack and fits_in_process:
            stack = fits_in_process[0]
    if stack is None:
        # Try root *.fit
        for p in session_dir.glob("*_og.fit"):
            stack = p
            break

    if stack and HAS_ASTROPY:
        with fits.open(stack) as hdul:
            h = hdul[0].header
            fi.ra = float(h.get("RA", h.get("CRVAL1", 0)))
            fi.dec = float(h.get("DEC", h.get("CRVAL2", 0)))
            fi.object_name = str(h.get("OBJECT", "Unknown"))
            fi.n_frames = int(h.get("STACKCNT", 0))
            fi.exptime = float(h.get("EXPTIME", 0))
            fi.date_obs = str(h.get("DATE-OBS", ""))
            fi.stack_path = stack
            try:
                fi.wcs = WCS(h)
            except Exception:
                fi.wcs = None

    # Load variable stars from starsv.csv if it exists
    stars: List[VariableStar] = []
    starsv_csv = session_dir / "starsv.csv"
    if starsv_csv.exists():
        with open(starsv_csv, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                try:
                    ra, dec = float(row["ra"]), float(row["dec"])
                    mag = float(row.get("mag") or row.get("magnitude") or 99)
                    name = row.get("name", row.get("Name", f"RA{ra:.3f}"))
                    # Distance from field center
                    dist = 0.0
                    if fi.ra and HAS_ASTROPY:
                        c_star = SkyCoord(ra, dec, unit=u.deg)
                        c_field = SkyCoord(fi.ra, fi.dec, unit=u.deg)
                        dist = c_star.separation(c_field).deg
                    stars.append(VariableStar(
                        name=name, ra=ra, dec=dec, mag=mag,
                        dist_from_center=round(dist, 3)
                    ))
                except (ValueError, KeyError):
                    continue

    return fi, stars


class VSXEnricher(QThread):
    """Queries VizieR VSX catalog to enrich variable star list with type/period."""
    progress = pyqtSignal(str)
    done = pyqtSignal(list)

    def __init__(self, stars: List[VariableStar], field: FieldInfo):
        super().__init__()
        self.stars = stars
        self.field = field

    def run(self):
        if not HAS_ASTROQUERY or not self.field.ra:
            self.done.emit(self.stars)
            return
        try:
            self.progress.emit("Querying AAVSO VSX catalog via VizieR...")
            coord = SkyCoord(self.field.ra, self.field.dec, unit=u.deg)
            # Search with generous radius to cover the full field
            radius = max(SEESTAR_PROFILE["fov_deg_w"], SEESTAR_PROFILE["fov_deg_h"])
            v = Vizier(columns=["Name", "Type", "Period", "max", "min", "u_max", "RAJ2000", "DEJ2000"],
                       row_limit=500)
            results = v.query_region(coord, radius=radius * u.deg, catalog="B/vsx/vsx")

            if not results or len(results) == 0:
                self.progress.emit("VSX: no results returned.")
                self.done.emit(self.stars)
                return

            vsx_table = results[0]
            self.progress.emit(f"VSX: {len(vsx_table)} entries in field.")

            # Build lookup dict by name (cleaned)
            def clean(s):
                return str(s).strip().lower().replace(" ", "").replace("*", "")

            vsx_by_name = {clean(r["Name"]): r for r in vsx_table}

            # Match and enrich
            for star in self.stars:
                key = clean(star.name)
                row = vsx_by_name.get(key)
                if row is None:
                    # Try coordinate match (within 30 arcsec)
                    for r in vsx_table:
                        try:
                            c = SkyCoord(float(r["RAJ2000"]), float(r["DEJ2000"]), unit=u.deg)
                            s = SkyCoord(star.ra, star.dec, unit=u.deg)
                            if c.separation(s).arcsec < 30:
                                row = r
                                break
                        except Exception:
                            continue
                if row is not None:
                    star.var_type = str(row.get("Type", "?")).strip() or "?"
                    per = row.get("Period")
                    star.period = f"{float(per):.3f}d" if per and str(per) != "--" else "—"
                    try:
                        mx = float(row["max"])
                        mn = float(row.get("min", mx))
                        star.amplitude = f"{abs(mn - mx):.2f} mag"
                    except Exception:
                        star.amplitude = "—"

            self.done.emit(self.stars)
        except Exception as exc:
            self.progress.emit(f"VSX query failed: {exc}")
            self.done.emit(self.stars)


class CompStarFetcher(QThread):
    """Queries APASS (VizieR II/336) for comparison stars near a target."""
    done = pyqtSignal(list)
    progress = pyqtSignal(str)

    def __init__(self, target: VariableStar, n_comp: int = 5):
        super().__init__()
        self.target = target
        self.n_comp = n_comp

    def run(self):
        if not HAS_ASTROQUERY:
            self.done.emit([])
            return
        try:
            self.progress.emit(f"Fetching comparison stars for {self.target.name} from APASS...")
            coord = SkyCoord(self.target.ra, self.target.dec, unit=u.deg)
            v = Vizier(
                columns=["RAJ2000", "DEJ2000", "Vmag", "Bmag", "e_Vmag"],
                column_filters={"e_Vmag": "<0.05"},
                row_limit=100,
            )
            result = v.query_region(coord, radius=1.0 * u.deg, catalog="II/336/apass9")
            if not result or len(result) == 0:
                self.progress.emit("APASS: no results. Using manual reference stars.")
                self.done.emit([])
                return

            tbl = result[0]
            comps: List[CompStar] = []
            target_vmag = self.target.mag

            for row in tbl:
                try:
                    vmag = float(row["Vmag"])
                    ra_c = float(row["RAJ2000"])
                    dec_c = float(row["DEJ2000"])
                    # Skip if too faint/bright vs target or at same position
                    if abs(vmag - target_vmag) > 2.0:
                        continue
                    c = SkyCoord(ra_c, dec_c, unit=u.deg)
                    t = SkyCoord(self.target.ra, self.target.dec, unit=u.deg)
                    if c.separation(t).arcsec < 5:
                        continue
                    comps.append(CompStar(ra=ra_c, dec=dec_c, vmag=vmag))
                except Exception:
                    continue

            # Sort by proximity to target magnitude, keep best N
            comps.sort(key=lambda c: abs(c.vmag - target_vmag))
            comps = comps[: self.n_comp]
            self.progress.emit(f"Found {len(comps)} comparison stars from APASS.")
            self.done.emit(comps)
        except Exception as exc:
            self.progress.emit(f"APASS query failed: {exc}")
            self.done.emit([])


class PipelineWorker(QThread):
    """Runs the full Siril preprocessing + photometry pipeline in background."""
    log = pyqtSignal(str)
    progress = pyqtSignal(int, str)   # percent, label
    done = pyqtSignal(bool, str)      # success, result_path

    def __init__(self, session_dir: Path, target: VariableStar,
                 comp_stars: List[CompStar], process_dir: Path):
        super().__init__()
        self.session_dir = session_dir
        self.target = target
        self.comp_stars = comp_stars
        self.process_dir = process_dir

    def _log(self, msg: str):
        self.log.emit(msg)

    def run(self):
        session = self.session_dir
        proc = self.process_dir
        lights = session / "lights"
        focal = SEESTAR_PROFILE["focal_mm"]
        pixsz = SEESTAR_PROFILE["pixel_um"]

        if not lights.is_dir():
            self.done.emit(False, f"lights/ folder not found in {session}")
            return

        n_fits = len([f for f in lights.iterdir() if f.suffix.lower() in (".fit", ".fits")])
        self._log(f"Found {n_fits} FITS frames in {lights}")

        # ── Step 1: Link lights into process/ ────────────────────────────────
        self.progress.emit(5, "Linking light frames...")
        self._log("── Step 1: link lights → sequence ──")
        cmds = [
            f'cd "{lights}"',
            f'link light -out="{proc}"',
            f'cd "{proc}"',
        ]
        if not run_siril(proc, cmds, self._log, "_step1_link.ssf"):
            self.done.emit(False, "Failed at step 1 (link)")
            return

        # link light → séquence nommée "light_" (sans 's')
        seq = "light_"

        # ── Step 2: Plate solve each frame ───────────────────────────────────
        self.progress.emit(20, "Plate solving sequence...")
        self._log("── Step 2: seqplatesolve ──")

        disto_dir = proc / "ps_distortion"
        disto_arg = f"-disto={disto_dir}" if disto_dir.is_dir() else ""

        cmds = [
            f'cd "{proc}"',
            # -nocache -force : recalcule les transformations (obligatoire pour seqapplyreg)
            f"seqplatesolve {seq} -nocache -force -focal={focal} -pixelsize={pixsz} "
            f"-radius=2.5 {disto_arg}".strip(),
        ]
        if not run_siril(proc, cmds, self._log, "_step2_platesolve.ssf"):
            self._log("WARNING: plate solve had errors — continuing anyway")

        # ── Step 3: Register frames ───────────────────────────────────────────
        self.progress.emit(55, "Registering frames...")
        self._log("── Step 3: seqapplyreg ──")
        cmds = [
            f'cd "{proc}"',
            f"seqapplyreg {seq} -framing=max -filter-round=2.5k",
        ]
        if not run_siril(proc, cmds, self._log, "_step3_reg.ssf"):
            self.done.emit(False, "Failed at step 3 (register)")
            return

        registered = f"r_{seq}"

        # ── Step 4: Photometry ────────────────────────────────────────────────
        self.progress.emit(80, "Running aperture photometry...")
        self._log("── Step 4: setphot + light_curve ──")

        lc_cmd = (
            f"light_curve {registered} {PHOT_CHANNEL} -autoring "
            f"-wcs={self.target.wcs_arg()}"
        )
        for cs in self.comp_stars:
            lc_cmd += f" -refwcs={cs.ra:.6f},{cs.dec:+.6f}"

        cmds = [
            f'cd "{proc}"',
            f"setphot -aperture={PHOT_APERTURE} -inner={PHOT_INNER} "
            f"-outer={PHOT_OUTER} -gain={SEESTAR_PROFILE['gain']}",
            lc_cmd,
        ]
        if not run_siril(proc, cmds, self._log, "_step4_phot.ssf"):
            self.done.emit(False, "Failed at step 4 (photometry)")
            return

        # ── Results ───────────────────────────────────────────────────────────
        lc_dat = proc / "light_curve.dat"
        results_dir = session / "results"
        results_dir.mkdir(exist_ok=True)

        safe_name = self.target.name.replace(" ", "_").replace("/", "-")
        out_dat = results_dir / f"{safe_name}_light_curve.dat"
        out_png = results_dir / f"{safe_name}_light_curve.png"
        out_csv = results_dir / f"{safe_name}_aavso.csv"

        if lc_dat.exists():
            shutil.copy2(lc_dat, out_dat)
            src_png = proc / "light_curve.png"
            if src_png.exists():
                shutil.copy2(src_png, out_png)
            _export_aavso(out_dat, out_csv, self.target.name)
            self.progress.emit(100, "Done!")
            self.done.emit(True, str(out_dat))
        else:
            self.done.emit(False, "light_curve.dat not found — check if target is in field")


def _export_aavso(lc_dat: Path, out_csv: Path, target_name: str):
    rows = []
    with open(lc_dat) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                vals = [float(p) for p in parts]
                # Detect JD column (> 2400000)
                jd_idx = next((i for i, v in enumerate(vals) if v > 2_400_000), 1)
                jd = vals[jd_idx]
                mag = vals[jd_idx + 1]
                err = vals[jd_idx + 2] if len(vals) > jd_idx + 2 else 0.0
                rows.append((jd, mag, err))
            except (ValueError, StopIteration):
                continue
    if not rows:
        return
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["#TYPE=EXTENDED"])
        w.writerow(["#OBSCODE=XXXX"])
        w.writerow(["#SOFTWARE=Siril + seestar_varstar_gui.py"])
        w.writerow(["#DELIM=,"])
        w.writerow(["#DATE=JD"])
        w.writerow(["#OBSTYPE=CCD"])
        w.writerow(["NAME", "DATE", "MAG", "MERR", "FILT", "TRANS", "MTYPE",
                    "CNAME", "CMAG", "KNAME", "KMAG", "AMASS", "GROUP", "CHART", "NOTES"])
        for jd, mag, err in rows:
            w.writerow([target_name, f"{jd:.6f}", f"{mag:.4f}", f"{err:.4f}",
                        "TG", "NO", "STD", "ENSEMBLE", "na", "na", "na",
                        "na", "1", "na", "seestar_s30pro"])


# ─────────────────────────────────────────────────────────────────────────────
# Table model
# ─────────────────────────────────────────────────────────────────────────────

COLS = ["Name", "Type", "Period", "Ampl.", "Mag", "Dist (°)"]

class StarTableModel(QAbstractTableModel):
    def __init__(self):
        super().__init__()
        self._stars: List[VariableStar] = []

    def set_stars(self, stars: List[VariableStar]):
        self.beginResetModel()
        self._stars = stars
        self.endResetModel()

    def star_at(self, row: int) -> Optional[VariableStar]:
        if 0 <= row < len(self._stars):
            return self._stars[row]
        return None

    def rowCount(self, parent=QModelIndex()):
        return len(self._stars)

    def columnCount(self, parent=QModelIndex()):
        return len(COLS)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            return COLS[section]
        return None

    def data(self, index: QModelIndex, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        star = self._stars[index.row()]
        col = index.column()
        if role == Qt.ItemDataRole.DisplayRole:
            return [star.name, star.var_type, star.period, star.amplitude,
                    f"{star.mag:.2f}", f"{star.dist_from_center:.2f}"][col]
        if role == Qt.ItemDataRole.ForegroundRole:
            if star.var_type not in ("?", ""):
                return QColor("#7EC8E3")   # known variable types in blue
        if role == Qt.ItemDataRole.TextAlignmentRole:
            if col in (2, 3, 4, 5):
                return Qt.AlignmentFlag.AlignCenter
        return None


# ─────────────────────────────────────────────────────────────────────────────
# Light curve plot widget
# ─────────────────────────────────────────────────────────────────────────────

class LightCurvePlot(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        if HAS_MATPLOTLIB:
            self.fig = Figure(figsize=(6, 3), facecolor="#1e1e2e")
            self.canvas = FigureCanvasQTAgg(self.fig)
            layout.addWidget(self.canvas)
            self.ax = self.fig.add_subplot(111)
            self._style_ax()
        else:
            layout.addWidget(QLabel("matplotlib not installed — no plot\npip3 install matplotlib"))

    def _style_ax(self):
        self.ax.set_facecolor("#13131f")
        self.ax.tick_params(colors="#aaaaaa")
        for spine in self.ax.spines.values():
            spine.set_color("#444466")
        self.ax.set_xlabel("JD", color="#aaaaaa", fontsize=8)
        self.ax.set_ylabel("Diff. magnitude", color="#aaaaaa", fontsize=8)
        self.ax.invert_yaxis()

    def load(self, dat_path: Path, star_name: str):
        if not HAS_MATPLOTLIB:
            return
        jds, mags, errs = [], [], []
        try:
            with open(dat_path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    vals = [float(p) for p in parts]
                    jd_idx = next((i for i, v in enumerate(vals) if v > 2_400_000), 1)
                    jds.append(vals[jd_idx])
                    mags.append(vals[jd_idx + 1])
                    errs.append(vals[jd_idx + 2] if len(vals) > jd_idx + 2 else 0.0)
        except Exception:
            return

        if not jds:
            return

        self.ax.clear()
        self._style_ax()
        jd0 = int(min(jds))
        x = [j - jd0 for j in jds]
        self.ax.errorbar(x, mags, yerr=errs, fmt="o", ms=3, color="#7EC8E3",
                         ecolor="#446688", capsize=2, lw=0.8)
        self.ax.set_title(f"{star_name}  (JD - {jd0})", color="#dddddd", fontsize=9)
        self.ax.set_xlabel(f"JD − {jd0}", color="#aaaaaa", fontsize=8)
        self.fig.tight_layout()
        self.canvas.draw()


# ─────────────────────────────────────────────────────────────────────────────
# Main window
# ─────────────────────────────────────────────────────────────────────────────

DARK_STYLE = """
QMainWindow, QWidget { background: #1e1e2e; color: #cdd6f4; font-size: 13px; }
QGroupBox { border: 1px solid #444466; border-radius: 4px; margin-top: 10px;
            padding-top: 8px; color: #89b4fa; font-weight: bold; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; }
QPushButton { background: #313244; border: 1px solid #45475a; border-radius: 4px;
              padding: 6px 14px; color: #cdd6f4; }
QPushButton:hover { background: #45475a; }
QPushButton#run_btn { background: #1e6a3e; border-color: #40a870; color: #ceffde;
                      font-weight: bold; font-size: 14px; padding: 8px 20px; }
QPushButton#run_btn:hover { background: #26844e; }
QPushButton#run_btn:disabled { background: #2a2a3a; color: #666; border-color: #444; }
QLineEdit { background: #181825; border: 1px solid #45475a; border-radius: 3px;
            padding: 4px 8px; color: #cdd6f4; }
QTableView { background: #181825; border: 1px solid #45475a; gridline-color: #2a2a3a;
             selection-background-color: #313244; alternate-background-color: #1a1a28; }
QTableView::item:selected { background: #45475a; color: #89b4fa; }
QHeaderView::section { background: #242436; color: #89b4fa; border: 1px solid #45475a;
                        padding: 4px 8px; font-weight: bold; }
QTextEdit { background: #11111b; border: 1px solid #45475a; border-radius: 3px;
            color: #a6adc8; font-family: "Menlo", "Courier New", monospace; font-size: 11px; }
QProgressBar { border: 1px solid #45475a; border-radius: 3px; background: #181825;
               text-align: center; color: #cdd6f4; }
QProgressBar::chunk { background: #1e6a3e; }
QLabel#field_label { color: #89dceb; font-size: 12px; }
QLabel#status_label { color: #a6e3a1; font-size: 12px; }
QComboBox { background: #313244; border: 1px solid #45475a; border-radius: 3px;
            padding: 4px 8px; color: #cdd6f4; }
"""


class MainWindow(QMainWindow):
    def __init__(self, initial_session: Optional[Path] = None):
        super().__init__()
        self.setWindowTitle("Seestar S30 Pro — Variable Star Finder & Photometry")
        self.setMinimumSize(1000, 720)
        self.setStyleSheet(DARK_STYLE)

        self.field_info: Optional[FieldInfo] = None
        self.all_stars: List[VariableStar] = []
        self.selected_star: Optional[VariableStar] = None
        self.comp_stars: List[CompStar] = []

        self._build_ui()
        self._check_siril()

        if initial_session:
            self.session_edit.setText(str(initial_session))
            self._load_session()
        else:
            self._auto_detect_session()

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(8)

        # ── Header ──
        header = QLabel("🔭  Seestar S30 Pro — Variable Star Finder & Photometry")
        header.setStyleSheet("font-size: 16px; font-weight: bold; color: #89b4fa; padding: 4px 0;")
        root.addWidget(header)

        # ── Session selector ──
        session_group = QGroupBox("Session")
        session_layout = QHBoxLayout(session_group)
        self.session_edit = QLineEdit()
        self.session_edit.setPlaceholderText("Path to session folder (e.g. /Volumes/Seestar/MyWorks/postprod_m81)")
        self.session_edit.returnPressed.connect(self._load_session)
        browse_btn = QPushButton("Browse…")
        browse_btn.clicked.connect(self._browse_session)
        load_btn = QPushButton("Load session")
        load_btn.clicked.connect(self._load_session)
        session_layout.addWidget(QLabel("Folder:"))
        session_layout.addWidget(self.session_edit, 1)
        session_layout.addWidget(browse_btn)
        session_layout.addWidget(load_btn)
        root.addWidget(session_group)

        # ── Field info bar ──
        field_group = QGroupBox("Field")
        field_layout = QHBoxLayout(field_group)
        self.lbl_object = QLabel("Object: —")
        self.lbl_ra = QLabel("RA: —")
        self.lbl_dec = QLabel("Dec: —")
        self.lbl_frames = QLabel("Frames: —")
        self.lbl_exptime = QLabel("Exp: —")
        self.lbl_siril = QLabel("")
        for lbl in [self.lbl_object, self.lbl_ra, self.lbl_dec,
                    self.lbl_frames, self.lbl_exptime]:
            lbl.setObjectName("field_label")
            field_layout.addWidget(lbl)
        field_layout.addStretch()
        field_layout.addWidget(self.lbl_siril)
        root.addWidget(field_group)

        # ── Main splitter: star table | details + log ──
        splitter = QSplitter(Qt.Orientation.Horizontal)
        root.addWidget(splitter, 1)

        # Left: star table
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)

        stars_group = QGroupBox("Variable stars in field")
        stars_layout = QVBoxLayout(stars_group)

        toolbar = QHBoxLayout()
        scan_btn = QPushButton("↻  Enrich with VSX catalog")
        scan_btn.clicked.connect(self._enrich_vsx)
        self.lbl_star_count = QLabel("0 stars")
        toolbar.addWidget(scan_btn)
        toolbar.addStretch()
        toolbar.addWidget(self.lbl_star_count)
        stars_layout.addLayout(toolbar)

        self.star_model = StarTableModel()
        self.proxy = QSortFilterProxyModel()
        self.proxy.setSourceModel(self.star_model)
        self.proxy.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        self.proxy.setFilterKeyColumn(-1)

        self.table = QTableView()
        self.table.setModel(self.proxy)
        self.table.setSortingEnabled(True)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableView.SelectionMode.SingleSelection)
        self.table.horizontalHeader().setStretchLastSection(False)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for i in range(1, len(COLS)):
            self.table.horizontalHeader().setSectionResizeMode(i, QHeaderView.ResizeMode.ResizeToContents)
        self.table.selectionModel().selectionChanged.connect(self._on_star_selected)
        stars_layout.addWidget(self.table)

        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Filter:"))
        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText("Search by name or type…")
        self.filter_edit.textChanged.connect(self.proxy.setFilterFixedString)
        filter_row.addWidget(self.filter_edit, 1)
        stars_layout.addLayout(filter_row)

        left_layout.addWidget(stars_group)
        splitter.addWidget(left)

        # Right: details + controls + log
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(8)

        # Target info
        target_group = QGroupBox("Selected target")
        target_layout = QGridLayout(target_group)
        self.lbl_target_name = QLabel("—")
        self.lbl_target_name.setStyleSheet("font-weight: bold; color: #f38ba8; font-size: 14px;")
        self.lbl_target_type = QLabel("—")
        self.lbl_target_coords = QLabel("—")
        self.lbl_target_period = QLabel("—")
        self.lbl_comp_count = QLabel("0 comparison stars")
        target_layout.addWidget(QLabel("Name:"), 0, 0)
        target_layout.addWidget(self.lbl_target_name, 0, 1)
        target_layout.addWidget(QLabel("Type / Period:"), 1, 0)
        target_layout.addWidget(self.lbl_target_type, 1, 1)
        target_layout.addWidget(QLabel("Coords:"), 2, 0)
        target_layout.addWidget(self.lbl_target_coords, 2, 1)
        target_layout.addWidget(QLabel("Comparison stars:"), 3, 0)
        target_layout.addWidget(self.lbl_comp_count, 3, 1)
        fetch_comp_btn = QPushButton("Fetch comparison stars (APASS)")
        fetch_comp_btn.clicked.connect(self._fetch_comp_stars)
        target_layout.addWidget(fetch_comp_btn, 4, 0, 1, 2)
        right_layout.addWidget(target_group)

        # Light curve plot
        plot_group = QGroupBox("Light curve")
        plot_layout = QVBoxLayout(plot_group)
        self.plot = LightCurvePlot()
        plot_layout.addWidget(self.plot)
        right_layout.addWidget(plot_group, 1)

        # Run button
        self.run_btn = QPushButton("▶  Prepare frames & Generate light curve")
        self.run_btn.setObjectName("run_btn")
        self.run_btn.setEnabled(False)
        self.run_btn.clicked.connect(self._run_pipeline)
        right_layout.addWidget(self.run_btn)

        # Progress
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_label = QLabel("Ready")
        self.progress_label.setObjectName("status_label")
        right_layout.addWidget(self.progress_bar)
        right_layout.addWidget(self.progress_label)

        # Log
        log_group = QGroupBox("Siril log")
        log_layout = QVBoxLayout(log_group)
        self.log_view = QTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumHeight(140)
        log_layout.addWidget(self.log_view)
        right_layout.addWidget(log_group)

        splitter.addWidget(right)
        splitter.setSizes([520, 480])

        # Status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("No session loaded")

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _check_siril(self):
        cli = find_siril_cli()
        if cli:
            self.lbl_siril.setText(f"✓ siril-cli: {cli}")
            self.lbl_siril.setStyleSheet("color: #a6e3a1;")
        else:
            self.lbl_siril.setText("✗ siril-cli not found")
            self.lbl_siril.setStyleSheet("color: #f38ba8;")

    def _log(self, msg: str):
        self.log_view.append(msg)
        self.log_view.verticalScrollBar().setValue(
            self.log_view.verticalScrollBar().maximum()
        )

    def _auto_detect_session(self):
        for base in SEESTAR_AUTO_PATHS:
            base_path = Path(base)
            if not base_path.exists():
                continue
            for child in sorted(base_path.iterdir(), reverse=True):
                if child.is_dir() and (child / "lights").is_dir():
                    self.session_edit.setText(str(child))
                    self._load_session()
                    return

    # ── Session loading ───────────────────────────────────────────────────────

    def _browse_session(self):
        d = QFileDialog.getExistingDirectory(self, "Select session folder",
                                             str(Path.home()))
        if d:
            self.session_edit.setText(d)
            self._load_session()

    def _load_session(self):
        path_str = self.session_edit.text().strip()
        if not path_str:
            return
        session = Path(path_str)
        if not session.is_dir():
            QMessageBox.warning(self, "Not found", f"Folder not found:\n{session}")
            return
        self._log(f"Loading session: {session}")
        self.status_bar.showMessage("Loading session…")
        self._loader = SessionLoader(session)
        self._loader.done.connect(self._on_session_loaded)
        self._loader.error.connect(lambda e: self._log(f"ERROR: {e}"))
        self._loader.start()

    def _on_session_loaded(self, fi: FieldInfo, stars: List[VariableStar]):
        self.field_info = fi
        self.all_stars = stars

        self.lbl_object.setText(f"Object: {fi.object_name}")
        self.lbl_ra.setText(f"RA: {fi.ra_hms}")
        self.lbl_dec.setText(f"Dec: {fi.dec_dms}")
        self.lbl_frames.setText(f"Stack: {fi.n_frames} × {fi.exptime:.0f}s")
        self.lbl_exptime.setText(f"Date: {fi.date_obs[:10] if fi.date_obs else '—'}")

        self.star_model.set_stars(stars)
        self.lbl_star_count.setText(f"{len(stars)} stars")
        self.status_bar.showMessage(
            f"Session loaded — {len(stars)} variable stars found in field"
        )
        self._log(f"Field: {fi.object_name}  RA={fi.ra:.4f}°  Dec={fi.dec:.4f}°")
        self._log(f"{len(stars)} variable stars loaded from starsv.csv")
        if not stars and fi.ra:
            self._log("Tip: click '↻ Enrich with VSX catalog' to query the online catalog.")

    # ── VSX enrichment ────────────────────────────────────────────────────────

    def _enrich_vsx(self):
        if not self.field_info:
            QMessageBox.information(self, "No session", "Load a session first.")
            return
        if not HAS_ASTROQUERY:
            QMessageBox.warning(self, "Missing dependency",
                                "astroquery not installed.\npip3 install astroquery")
            return
        self.status_bar.showMessage("Querying VSX catalog…")
        self._vsx_worker = VSXEnricher(self.all_stars, self.field_info)
        self._vsx_worker.progress.connect(self._log)
        self._vsx_worker.done.connect(self._on_vsx_done)
        self._vsx_worker.start()

    def _on_vsx_done(self, stars: List[VariableStar]):
        self.all_stars = stars
        self.star_model.set_stars(stars)
        self.status_bar.showMessage(f"VSX enrichment done — {len(stars)} stars")
        self._log("VSX enrichment complete.")

    # ── Star selection ────────────────────────────────────────────────────────

    def _on_star_selected(self):
        indices = self.table.selectionModel().selectedRows()
        if not indices:
            self.selected_star = None
            self.run_btn.setEnabled(False)
            return
        src_idx = self.proxy.mapToSource(indices[0])
        star = self.star_model.star_at(src_idx.row())
        if not star:
            return
        self.selected_star = star
        self.comp_stars = []

        self.lbl_target_name.setText(star.name)
        self.lbl_target_type.setText(
            f"{star.var_type}  ·  Period: {star.period}  ·  Ampl.: {star.amplitude}"
        )
        self.lbl_target_coords.setText(
            f"RA {star.ra:.6f}°  /  Dec {star.dec:+.6f}°"
        )
        self.lbl_comp_count.setText("0 comparison stars (click 'Fetch' to load)")
        self.run_btn.setEnabled(True)
        self.status_bar.showMessage(f"Selected: {star.name}")

    # ── Comparison stars ──────────────────────────────────────────────────────

    def _fetch_comp_stars(self):
        if not self.selected_star:
            return
        if not HAS_ASTROQUERY:
            QMessageBox.warning(self, "Missing", "pip3 install astroquery")
            return
        self._comp_worker = CompStarFetcher(self.selected_star)
        self._comp_worker.progress.connect(self._log)
        self._comp_worker.done.connect(self._on_comp_done)
        self._comp_worker.start()

    def _on_comp_done(self, comps: List[CompStar]):
        self.comp_stars = comps
        self.lbl_comp_count.setText(f"{len(comps)} comparison stars from APASS")
        self._log(f"Comparison stars ready: {len(comps)}")

    # ── Pipeline ──────────────────────────────────────────────────────────────

    def _run_pipeline(self):
        if not self.selected_star:
            return
        if not self.field_info:
            return

        session = Path(self.session_edit.text().strip())
        proc_dir = session / "process"
        proc_dir.mkdir(exist_ok=True)

        if len(self.comp_stars) < 3:
            ans = QMessageBox.question(
                self, "Few comparison stars",
                f"Only {len(self.comp_stars)} comparison stars loaded (minimum recommended: 3).\n"
                "Click 'Fetch comparison stars (APASS)' first for better accuracy.\n\n"
                "Continue anyway?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            if ans != QMessageBox.StandardButton.Yes:
                return

        self.run_btn.setEnabled(False)
        self.progress_bar.setValue(0)
        self.log_view.clear()
        self._log(f"Starting pipeline for: {self.selected_star.name}")
        self._log(f"Target: RA={self.selected_star.ra:.6f}  Dec={self.selected_star.dec:.6f}")
        self._log(f"Comparison stars: {len(self.comp_stars)}")

        self._pipeline = PipelineWorker(session, self.selected_star, self.comp_stars, proc_dir)
        self._pipeline.log.connect(self._log)
        self._pipeline.progress.connect(self._on_progress)
        self._pipeline.done.connect(self._on_pipeline_done)
        self._pipeline.start()

    def _on_progress(self, pct: int, label: str):
        self.progress_bar.setValue(pct)
        self.progress_label.setText(label)
        self.status_bar.showMessage(label)

    def _on_pipeline_done(self, success: bool, result: str):
        self.run_btn.setEnabled(True)
        if success:
            result_path = Path(result)
            self._log(f"✓ Light curve saved: {result_path}")
            self.status_bar.showMessage(f"Done! Results in {result_path.parent}")
            QMessageBox.information(self, "Done",
                                    f"Light curve generated!\n\n{result_path}\n\n"
                                    f"AAVSO CSV: {result_path.parent / (result_path.stem.replace('_light_curve', '_aavso') + '.csv')}")
            if result_path.exists():
                self.plot.load(result_path, self.selected_star.name)
        else:
            self._log(f"✗ Pipeline failed: {result}")
            self.status_bar.showMessage(f"Failed: {result}")
            QMessageBox.critical(self, "Pipeline failed", result)


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Seestar variable star photometry GUI")
    parser.add_argument("--session", help="Session folder path (auto-detected if omitted)")
    args = parser.parse_args()

    app = QApplication(sys.argv)
    app.setApplicationName("Seestar VarStar")

    session_path = Path(args.session) if args.session else None
    win = MainWindow(initial_session=session_path)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
