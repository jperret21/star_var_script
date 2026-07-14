"""
pipeline.py — Seestar S30 Pro variable-star photometry pipeline.

Pure logic: FITS I/O, catalog queries, Siril orchestration, AAVSO export.
No tkinter dependency — importable from tests, CLI scripts, or the GUI.
"""

from __future__ import annotations

import csv
import datetime
import math
import os
import platform
import shutil
import statistics
import subprocess
from pathlib import Path
from typing import Callable, Optional

VERSION = "0.1.5"

# ── optional runtime deps ────────────────────────────────────────────────────

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

try:
    import sys as _sys
    _sys.path.insert(0, "/Applications/Siril.app/Contents/Resources/share/siril/python_module")
    import sirilpy
    from sirilpy import SirilInterface
    HAS_SIRILPY = True
except ImportError:
    HAS_SIRILPY = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# ── hardware profile ─────────────────────────────────────────────────────────

SEESTAR = {
    "focal":  160,    # mm  (FOCALLEN from FITS)
    "pixsz":  2.9,    # µm  (XPIXSZ from FITS)
    "fov_w":  2.2,    # deg width  (~2160px × 2.9µm / 160mm)
    "fov_h":  3.9,    # deg height (~3840px × 2.9µm / 160mm)
    # gain: read from FITS EGAIN header at runtime (GAIN=200 is ISO, not e-/ADU)
}

VIZIER = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"

# AAVSO filter codes — (code, extra_note)
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
# FITS utilities (pure Python, no astropy)
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


def _fits_strip_keyword(path: Path, keyword: str) -> Optional[tuple[int, str]]:
    """Remove a FITS header keyword in-place. Returns (byte_offset, value) for restore."""
    kw_bytes = keyword.upper().ljust(8, ' ').encode('ascii')
    try:
        with open(path, 'r+b') as f:
            block_start = 0
            while True:
                block = f.read(2880)
                if len(block) < 80:
                    break
                for i in range(0, len(block), 80):
                    card = block[i:i+80]
                    if card[:8] == kw_bytes:
                        val = card[10:80].decode('ascii', errors='ignore')
                        val = val.split('/')[0].strip().strip("'").strip()
                        f.seek(block_start + i)
                        f.write(b' ' * 80)
                        return (block_start + i, val)
                    if card[:3] == b'END' and (len(card) < 4 or card[3:4] in (b' ', b'\x00')):
                        return None
                block_start += 2880
    except Exception:
        pass
    return None


def _fits_restore_keyword(path: Path, keyword: str, offset: int, value: str) -> bool:
    """Restore a FITS keyword card at the exact byte offset it was removed from."""
    kw = keyword.upper().ljust(8, ' ')[:8]
    card = f"{kw}= '{value:<8}'"
    card_bytes = card.encode('ascii').ljust(80)[:80]
    try:
        with open(path, 'r+b') as f:
            f.seek(offset)
            f.write(card_bytes)
        return True
    except Exception:
        return False


def sky_to_pixel(ra: float, dec: float, hdr: dict):
    """
    TAN gnomonic projection (ignores SIP, which is sub-pixel at Seestar scale).
    Returns (px, py) in 1-based FITS pixel coords, or (None, None) on failure.
    """
    try:
        crval1 = float(hdr["CRVAL1"]); crval2 = float(hdr["CRVAL2"])
        crpix1 = float(hdr["CRPIX1"]); crpix2 = float(hdr["CRPIX2"])
        if "CD1_1" in hdr:
            cd11 = float(hdr["CD1_1"]); cd12 = float(hdr["CD1_2"])
            cd21 = float(hdr["CD2_1"]); cd22 = float(hdr["CD2_2"])
        else:
            cdelt1 = float(hdr["CDELT1"]); cdelt2 = float(hdr["CDELT2"])
            cd11 = cdelt1 * float(hdr.get("PC1_1", 1.0))
            cd12 = cdelt1 * float(hdr.get("PC1_2", 0.0))
            cd21 = cdelt2 * float(hdr.get("PC2_1", 0.0))
            cd22 = cdelt2 * float(hdr.get("PC2_2", 1.0))
        ra0 = math.radians(crval1); dec0 = math.radians(crval2)
        ra_r = math.radians(ra);    dec_r = math.radians(dec)
        dra = ra_r - ra0
        denom = (math.sin(dec0) * math.sin(dec_r) +
                 math.cos(dec0) * math.cos(dec_r) * math.cos(dra))
        if abs(denom) < 1e-10:
            return None, None
        x = math.degrees(math.cos(dec_r) * math.sin(dra) / denom)
        y = math.degrees((math.cos(dec0) * math.sin(dec_r) -
                          math.sin(dec0) * math.cos(dec_r) * math.cos(dra)) / denom)
        det = cd11 * cd22 - cd12 * cd21
        if abs(det) < 1e-20:
            return None, None
        px = crpix1 + (cd22 * x - cd12 * y) / det
        py = crpix2 + (-cd21 * x + cd11 * y) / det
        return px, py
    except (KeyError, ValueError, ZeroDivisionError):
        return None, None


def stars_in_frame(stars: list[dict], wcs_hdr: dict,
                   naxis1: int, naxis2: int,
                   margin: int = 50) -> list[dict]:
    """Return only stars whose sky coords project inside the image frame.

    margin: minimum pixel distance from any edge. Set to at least the outer
    photometry ring (30 px) so apertures don't clip at the boundary.
    """
    inside = []
    for star in stars:
        px, py = sky_to_pixel(star["ra"], star["dec"], wcs_hdr)
        if px is None:
            continue
        if (margin < px <= naxis1 - margin and
                margin < py <= naxis2 - margin):
            inside.append(star)
    return inside


def stars_in_safe_circle(stars: list[dict], wcs_hdr: dict,
                          naxis1: int, naxis2: int,
                          margin: int = 50) -> list[dict]:
    """Return stars inside the inscribed circle of the frame.

    For alt-az mounts, field rotation means only the inscribed circle
    (radius = min(naxis1, naxis2) / 2) is guaranteed to be covered by
    every registered frame.  Stars outside this circle may land in the
    black rotation corners of some frames.

    margin: pixels of safety buffer inside the circle boundary (absorbs
    registration drift and photometry aperture clipping).
    """
    if not wcs_hdr:
        return stars
    try:
        crpix1 = float(wcs_hdr.get("CRPIX1", (naxis1 + 1) / 2.0))
        crpix2 = float(wcs_hdr.get("CRPIX2", (naxis2 + 1) / 2.0))
    except (ValueError, TypeError):
        return stars
    safe_r = min(naxis1, naxis2) / 2.0 - margin
    if safe_r <= 0:
        return []
    safe_r2 = safe_r ** 2
    inside = []
    for star in stars:
        px, py = sky_to_pixel(star["ra"], star["dec"], wcs_hdr)
        if px is None:
            continue
        if (px - crpix1) ** 2 + (py - crpix2) ** 2 <= safe_r2:
            inside.append(star)
    return inside


def find_siril_user_catalogue() -> Optional[Path]:
    """Return path to Siril's user-DSO-catalogue.csv, or None if not found."""
    system = platform.system()
    if system == "Darwin":
        candidates = [
            Path.home() / "Library/Application Support/org.siril.Siril/siril/catalogue/user-DSO-catalogue.csv",
        ]
    elif system == "Linux":
        candidates = [
            Path.home() / ".local/share/siril/catalogue/user-DSO-catalogue.csv",
            Path.home() / ".siril/catalogue/user-DSO-catalogue.csv",
        ]
    else:  # Windows
        appdata = Path(os.environ.get("APPDATA", Path.home()))
        candidates = [appdata / "siril/catalogue/user-DSO-catalogue.csv"]

    for p in candidates:
        if p.parent.is_dir():
            return p
    return None


def update_siril_catalogue(stars: list[dict], catalogue_path: Path) -> int:
    """Append in-frame VSX stars to Siril's user DSO catalogue.

    Skips stars already present (matched by name or alias, both with and
    without the 'V* ' prefix that VSX sometimes adds).
    Returns the number of new entries written.
    """
    existing: set[str] = set()
    if catalogue_path.exists():
        with open(catalogue_path, encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                for field in ("name", "alias"):
                    v = row.get(field, "").strip()
                    if v:
                        existing.add(v)
                        existing.add(v.removeprefix("V* "))
                        existing.add(f"V* {v}")

    new_stars = [
        s for s in stars
        if s["name"] not in existing
        and s["name"].removeprefix("V* ") not in existing
    ]
    if not new_stars:
        return 0

    write_header = not catalogue_path.exists() or catalogue_path.stat().st_size == 0
    with open(catalogue_path, "a", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        if write_header:
            w.writerow(["name", "ra", "dec", "pmra", "pmdec", "mag", "bmag", "alias"])
        for star in new_stars:
            w.writerow([
                star["name"],
                f"{star['ra']:.6f}",
                f"{star['dec']:.6f}",
                "",
                "",
                f"{star.get('mag', 99.0):.2f}",
                "",
                star.get("var_type", ""),
            ])
    return len(new_stars)


# ─────────────────────────────────────────────────────────────────────────────
# Catalog queries via VizieR HTTP (no astroquery)
# ─────────────────────────────────────────────────────────────────────────────

def _vizier_tsv(catalog: str, ra: float, dec: float,
                radius_arcmin: float, columns: list[str],
                filters: dict | None = None, max_rows: int = 500) -> list[list[str]]:
    """Query VizieR and return parsed data rows. Returns [] on any error."""
    if not HAS_REQUESTS:
        return []
    params: dict = {
        "-source": catalog,
        "-c": f"{ra:.6f} {dec:+.6f}",
        "-c.r": str(radius_arcmin),
        "-c.u": "arcmin",
        "-out": ",".join(columns),
        "-out.max": str(max_rows),
        # Sort by distance from field centre: without this, VizieR returns rows
        # in catalog order and -out.max can truncate away the central stars
        # (e.g. RR Lyr lost among >500 KIC variables in the Kepler field).
        "-sort": "_r",
    }
    if filters:
        params.update(filters)
    try:
        r = requests.get(VIZIER, params=params, timeout=20)
        r.raise_for_status()
        # TSV layout: #comments, column names, units, dashed separator, data.
        rows = []
        in_data = False
        for line in r.text.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith("-"):
                in_data = True
                continue
            if in_data:
                rows.append(line.split("\t"))
        return rows
    except Exception:
        return []


def query_vsx(ra: float, dec: float, radius_deg: float) -> list[dict]:
    """Variable stars from VizieR B/vsx/vsx."""
    rows = _vizier_tsv(
        "B/vsx/vsx", ra, dec, radius_deg * 60,
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
                radius_deg: float = 1.0, n: int = 6,
                dvmag: float = 2.0) -> list[dict]:
    """Comparison stars from VizieR II/336/apass9.

    dvmag: max magnitude difference from target. Use 3.5+ when target_mag
    comes from VSX (which stores max brightness, often 1-2 mag brighter than mean).
    """
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
            if abs(vmag - target_mag) > dvmag:
                continue
            sep = ((ra_c - ra) * 3600) ** 2 + ((dec_c - dec) * 3600) ** 2
            if sep < 25:  # < 5 arcsec from target
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
    def __init__(self, log_cb: Callable[[str], None] | None = None):
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


# ─────────────────────────────────────────────────────────────────────────────
# Pipeline utilities
# ─────────────────────────────────────────────────────────────────────────────

def truncate_comp_csv(csv_path: Path, n: int) -> int:
    """Keep only the first n Comp1 entries in a findcompstars CSV."""
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


def _parse_comp_csv_vmag(path: Path) -> Optional[float]:
    """Return median V mag of comparison stars from a Siril findcompstars CSV.

    Format: comment lines (#), then 'type,name,ra,dec,mag' header, then data rows
    where rows starting with 'Comp' have catalog V magnitudes in column 5 (index 4).
    """
    vmags: list[float] = []
    try:
        with open(path, newline="", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(",")
                if len(parts) >= 5 and parts[0].startswith("Comp"):
                    try:
                        vmags.append(float(parts[4]))
                    except ValueError:
                        pass
    except Exception:
        return None
    return statistics.median(vmags) if vmags else None


def export_fwhm_csv(seq_path: Path, proc: Path, stem: str,
                    fixlen: int, out: Path,
                    plate_scale_arcsec: float = 1.035) -> int:
    """Extract per-frame FWHM from a Siril .seq file and write fwhm.csv.

    The .seq R0 lines contain FWHM in pixels from the registration star detection.
    Each R0 line is matched to its frame via the corresponding I line.
    DATE-OBS is read from each FITS header to build the JD axis.

    Returns number of rows written (0 on failure).
    """
    # Parse .seq: build list of (frame_num, selected, fwhm_x, fwhm_y)
    entries: list[tuple[int, bool, float, float]] = []
    try:
        frame_nums: list[int] = []
        selected:   list[bool] = []
        with open(seq_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("I "):
                    parts = line.split()
                    frame_nums.append(int(parts[1]))
                    selected.append(parts[2].strip() == "1")
        # R0 lines appear in the same order as I lines
        r0_idx = 0
        r0_data: list[tuple[float, float]] = []
        with open(seq_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("R0 "):
                    parts = line.split()
                    # R0 fwhm_x fwhm_y roundness ...
                    r0_data.append((float(parts[1]), float(parts[2])))
        for i, (fn, sel) in enumerate(zip(frame_nums, selected)):
            if i < len(r0_data):
                fx, fy = r0_data[i]
                entries.append((fn, sel, fx, fy))
    except Exception:
        return 0

    if not entries:
        return 0

    # Read DATE-OBS from each selected FITS header → JD
    rows: list[tuple[float, float, float]] = []   # (jd, fwhm_x_as, fwhm_y_as)
    for fn, sel, fx, fy in entries:
        if not sel or fx <= 0:
            continue
        fit = proc / f"{stem}_{fn:0{fixlen}d}.fit"
        if not fit.exists():
            continue
        hdr  = read_fits_header(fit)
        dobs = hdr.get("DATE-OBS") or hdr.get("DATE_OBS") or hdr.get("DATE")
        if not dobs:
            continue
        jd = dateobs_to_jd(str(dobs))
        if jd is None:
            continue
        rows.append((jd, fx * plate_scale_arcsec, fy * plate_scale_arcsec))

    if not rows:
        return 0

    rows.sort(key=lambda r: r[0])
    try:
        with open(out, "w", newline="", encoding="utf-8") as f:
            f.write("# FWHM per frame from Siril registration star detection\n")
            f.write(f"# Plate scale: {plate_scale_arcsec:.4f} arcsec/px\n")
            f.write("JD,FWHM_x_arcsec,FWHM_y_arcsec\n")
            for jd, fx, fy in rows:
                f.write(f"{jd:.6f},{fx:.3f},{fy:.3f}\n")
    except Exception:
        return 0
    return len(rows)


def dateobs_to_jd(date_str: str) -> Optional[float]:
    """Convert a FITS DATE-OBS string (ISO 8601) to Julian Date."""
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
    """Normalize light_curve.dat to absolute JD regardless of Siril's output format.

    Siril 1.4.3 outputs three possible formats depending on whether date_obs was
    populated in the sequence:
      - '#JD_UT (+ 0)' + frame indices 1,2,3…  → inject JD from FITS DATE-OBS
      - '#JD_UT (+ N)' + fractional offsets     → add N to convert to absolute JD
      - '#JD_UT'        + absolute JD already   → already correct, no-op

    Returns True if file was modified.
    """
    try:
        lines = dat_path.read_text(encoding="utf-8").splitlines()
    except Exception:
        return False

    julian0 = None
    for line in lines:
        if line.startswith("#JD_UT"):
            if "(+" in line:
                try:
                    julian0 = int(line.split("(+")[1].split(")")[0].strip())
                except (ValueError, IndexError):
                    julian0 = 0
            else:
                julian0 = None  # already absolute JD
            break

    if julian0 is None:
        return False

    out: list[str] = []

    if julian0 > 2_400_000:
        # Correct Siril output: fractional day offsets from julian0
        for line in lines:
            if line.startswith("#JD_UT"):
                out.append("#JD_UT")
                continue
            if line.startswith("#"):
                out.append(line)
                continue
            parts = line.split()
            if parts:
                try:
                    offset = float(parts[0])
                    if offset < 10:
                        parts[0] = f"{julian0 + offset:.6f}"
                except (ValueError, IndexError):
                    pass
            out.append(" ".join(parts))
        dat_path.write_text("\n".join(out) + "\n", encoding="utf-8")
        return True

    # julian0 == 0: frame indices → inject from FITS DATE-OBS
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

    jd_map: dict[int, float] = {}
    half_exp = exptime_s / 86400.0 / 2.0
    for idx, img_num in enumerate(selected_imgs, start=1):
        fit_path = proc / f"{stem}_{img_num:0{fixlen}d}.fit"
        hdr = read_fits_header(fit_path)
        date_str = str(hdr.get("DATE-OBS", ""))
        if date_str:
            jd = dateobs_to_jd(date_str)
            if jd is not None:
                jd_map[idx] = jd + half_exp

    if not jd_map:
        return False

    for line in lines:
        if line.startswith("#JD_UT"):
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


# ─────────────────────────────────────────────────────────────────────────────
# Frame-directory resolution (tolerant to capture-tool layouts)
# ─────────────────────────────────────────────────────────────────────────────

_FITS_GLOBS = ("*.fit", "*.fits", "*.fts")


def _dir_has_fits(d: Path) -> bool:
    """True if *d* directly contains at least one FITS file."""
    return d.is_dir() and any(next(iter(d.glob(g)), None) for g in _FITS_GLOBS)


def list_fits(d: Path) -> list[Path]:
    """All FITS files directly inside *d* (any common extension), sorted by name."""
    out: list[Path] = []
    for g in _FITS_GLOBS:
        out += d.glob(g)
    return sorted(out)


def resolve_frame_dir(session: Path, names) -> Optional[Path]:
    """Locate the directory holding frames of a given type under *session*.

    Tolerant to two capture layouts:
      • flat   — session/<name>/*.fits              (legacy Seestar dump)
      • nested — session/<Name>/<subfolder>/*.fits  (Argos: a per-filter
                 subfolder, e.g. Lights/IR/, Darks/Dark/, Biases/Dark/)

    *names* is a single folder name or a list of accepted names, matched
    case-insensitively (e.g. ["lights", "light"]). When the top folder holds
    frames directly it is returned as-is; otherwise the populated subfolder is
    returned — the richest one when several exist (multiple filters). Returns
    None when nothing is found.
    """
    if isinstance(names, str):
        names = [names]
    wanted = {n.lower() for n in names}
    if not session.is_dir():
        return None

    top = next(
        (c for c in sorted(session.iterdir())
         if c.is_dir() and c.name.lower() in wanted),
        None,
    )
    if top is None:
        return None
    if _dir_has_fits(top):
        return top

    subs = [s for s in sorted(top.iterdir()) if _dir_has_fits(s)]
    if not subs:
        return None
    return max(subs, key=lambda s: len(list_fits(s)))


def fits_exptime(d: Path) -> Optional[float]:
    """Exposure time (s) of the first FITS in *d*, or None if unreadable."""
    for fit in list_fits(d)[:1]:
        hdr = read_fits_header(fit)
        for key in ("EXPTIME", "EXPOSURE"):
            v = hdr.get(key)
            if v is not None:
                try:
                    return float(v)
                except (ValueError, TypeError):
                    pass
    return None


def get_gain_eadu(lights_dir: Path) -> float:
    """Read real gain in e-/ADU from the first FITS file in lights/.
    GAIN=200 on a Seestar is ISO-equivalent, not e-/ADU.
    Returns 1.0 as a safe fallback.
    """
    for fit in list_fits(lights_dir)[:1]:
        hdr = read_fits_header(fit)
        for key in ("EGAIN", "EPERDN", "GAIN_E", "CCDGAIN"):
            v = hdr.get(key)
            if v is not None:
                try:
                    g = float(v)
                    if 0.05 < g < 30:
                        return round(g, 3)
                except (ValueError, TypeError):
                    pass
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
# AAVSO export
# ─────────────────────────────────────────────────────────────────────────────

def build_airmass_map(proc: Path, stem: str, fixlen: int) -> list[tuple[float, float]]:
    """Build a per-frame (JD, airmass) table from the registered FITS headers.

    The Seestar/Argos capture tool writes an AIRMASS keyword in each frame header;
    DATE-OBS gives the JD axis. Returns a list sorted by JD, used to fill the AAVSO
    AMASS column by nearest-JD match. Empty if no frame exposes AIRMASS.
    """
    amap: list[tuple[float, float]] = []
    for fit in sorted(proc.glob(f"{stem}_*.fit")):
        hdr = read_fits_header(fit)
        am  = hdr.get("AIRMASS")
        if am is None:
            continue
        try:
            am = float(am)
        except (ValueError, TypeError):
            continue
        dobs = hdr.get("DATE-OBS") or hdr.get("DATE_OBS")
        if not dobs:
            continue
        jd = dateobs_to_jd(str(dobs))
        if jd is not None:
            amap.append((jd, am))
    amap.sort()
    return amap


def _airmass_at(jd: float, amap: list[tuple[float, float]],
                tol_days: float) -> Optional[float]:
    """Nearest airmass to *jd* within *tol_days*, or None if no frame is close enough."""
    if not amap:
        return None
    best_jd, best_am = min(amap, key=lambda r: abs(r[0] - jd))
    return best_am if abs(best_jd - jd) <= tol_days else None


def _robust_sigma(values: list[float]) -> float:
    """MAD-based robust standard deviation (1.4826 × median absolute deviation)."""
    if not values:
        return 0.0
    med = statistics.median(values)
    mad = statistics.median([abs(v - med) for v in values])
    return 1.4826 * mad


def estimate_lightcurve_scatter(jd_mag: list[tuple[float, float]]) -> Optional[float]:
    """Estimate the real per-point photometric scatter of a light curve, in mag.

    Uses robust second differences d2 = m[i-1] - 2·m[i] + m[i+1], which cancel any
    locally-linear trend, so a smoothly-varying star (e.g. RR Lyr) doesn't bias the
    estimate. Independent per-point noise σ propagates as Var(d2)=6σ², hence σ =
    robust_sigma(d2)/√6. Residual curvature only *inflates* the estimate, so the
    result is conservative. Returns σ (mag), or None if fewer than 10 usable points.
    """
    pts = sorted((jd, m) for jd, m in jd_mag
                 if not math.isnan(jd) and not math.isnan(m))
    if len(pts) < 10:
        return None
    d2 = [pts[i - 1][1] - 2 * pts[i][1] + pts[i + 1][1]
          for i in range(1, len(pts) - 1)]
    return _robust_sigma(d2) / math.sqrt(6.0)


def apply_error_floor(dat: Path,
                      manual_floor: Optional[float] = None) -> Optional[dict]:
    """Lift the error column of a light_curve.dat to a realistic value, in place.

    Siril reports only the formal photon-noise error, which underestimates the true
    scatter (flat-field residuals, scintillation, imperfect dark/flat calibration —
    none of which enter the photon budget). This measures the actual scatter σ_real
    and rewrites each point's error as (model B, quadrature):

        err' = √(err² + σ_syst²),   σ_syst = √(max(0, σ_real² − median(err)²))

    so every point keeps its own photon noise and gains a shared systematic floor.

    If *manual_floor* is given it is used directly as σ_syst (measurement skipped).
    Returns {sigma_real, median_err, sigma_syst, n} or None if not applied.
    """
    try:
        lines = dat.read_text(encoding="utf-8").splitlines()
    except Exception:
        return None

    # (line_index, err_col, jd, mag, err, parts) for each data row
    parsed: list[tuple[int, int, float, float, float, list]] = []
    for idx, line in enumerate(lines):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        try:
            jd_i = next(i for i, p in enumerate(parts) if float(p) > 2_400_000)
            jd   = float(parts[jd_i])
            mag  = float(parts[jd_i + 1])
            err  = float(parts[jd_i + 2])
        except (ValueError, StopIteration, IndexError):
            continue
        parsed.append((idx, jd_i + 2, jd, mag, err, parts))

    good = [(jd, mag, err) for (_i, _e, jd, mag, err, _p) in parsed
            if not math.isnan(mag) and not math.isnan(err)]
    if not good or (len(good) < 10 and manual_floor is None):
        return None

    med_err = statistics.median([e for _j, _m, e in good])
    if manual_floor is not None:
        sigma_syst = max(0.0, float(manual_floor))
        sigma_real = math.sqrt(med_err ** 2 + sigma_syst ** 2)
    else:
        sigma_real = estimate_lightcurve_scatter([(j, m) for j, m, _e in good])
        if sigma_real is None:
            return None
        sigma_syst = math.sqrt(max(0.0, sigma_real ** 2 - med_err ** 2))

    if sigma_syst > 0:
        for idx, err_col, jd, mag, err, parts in parsed:
            if math.isnan(err) or err_col >= len(parts):
                continue
            parts[err_col] = f"{math.sqrt(err ** 2 + sigma_syst ** 2):.6g}"
            lines[idx] = " ".join(parts)
        note = (f"# MERR includes systematic error floor: sigma_syst={sigma_syst:.4f} mag "
                f"(measured scatter={sigma_real:.4f}, formal median={med_err:.4f}; "
                f"err = sqrt(formal^2 + sigma_syst^2))")
        lines.insert(1 if lines and lines[0].startswith("#") else 0, note)
        dat.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return {"sigma_real": sigma_real, "median_err": med_err,
            "sigma_syst": sigma_syst, "n": len(good)}


def export_aavso(dat: Path, out: Path, name: str,
                 filt_code: str = "CV", filt_note: str = "",
                 merr_max: float = 0.5,
                 airmass_map: Optional[list[tuple[float, float]]] = None) -> list:
    """Export light curve to AAVSO Extended format.

    Filters NaN magnitudes and MERR > merr_max.
    When *airmass_map* (a sorted list of (JD, airmass) from build_airmass_map) is
    provided, the AMASS column is filled by nearest-JD match; otherwise it is "na".
    Returns list of (jd, mag, err) rows written.
    """
    rows = []
    with open(dat) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                jd_i = next(i for i, p in enumerate(parts) if float(p) > 2_400_000)
                jd   = float(parts[jd_i])
                mag  = float(parts[jd_i + 1])
                err  = float(parts[jd_i + 2]) if len(parts) > jd_i + 2 else 0.0
                if math.isnan(mag) or math.isnan(err):
                    continue
                if err > merr_max:
                    continue
                rows.append((jd, mag, err))
            except (ValueError, StopIteration, IndexError):
                continue
    if not rows:
        return []

    # Airmass match tolerance: half the median frame cadence, so each light-curve
    # point maps to exactly one frame. Falls back to 60 s for a sparse/1-point map.
    airmass_tol = 60.0 / 86400.0
    if airmass_map and len(airmass_map) >= 2:
        spacings = [airmass_map[i + 1][0] - airmass_map[i][0]
                    for i in range(len(airmass_map) - 1)]
        airmass_tol = statistics.median(spacings) / 2.0

    notes = f"seestar_s30pro|{filt_note}" if filt_note else "seestar_s30pro"
    with open(out, "w", newline="") as f:
        # Write header directly — csv.writer would quote "#DELIM=," (contains comma)
        for line in ["#TYPE=EXTENDED", "#OBSCODE=XXXX",
                     f"#SOFTWARE=Siril+seestar_varstar_siril.py v{VERSION}",
                     f"#FILTER={filt_code}",
                     "#DELIM=,", "#DATE=JD", "#OBSTYPE=CCD"]:
            f.write(line + "\n")
        w = csv.writer(f)
        w.writerow(["NAME","DATE","MAG","MERR","FILT","TRANS","MTYPE",
                    "CNAME","CMAG","KNAME","KMAG","AMASS","GROUP","CHART","NOTES"])
        for jd, mag, err in rows:
            am = _airmass_at(jd, airmass_map, airmass_tol) if airmass_map else None
            amass = f"{am:.4f}" if am is not None else "na"
            w.writerow([name, f"{jd:.6f}", f"{mag:.4f}", f"{err:.4f}",
                        filt_code, "NO", "DIFF",
                        "ENSEMBLE", "na", "na", "na", amass, "1", "na", notes])
    return rows


def export_photometry_csv(dat: Path, out: Path, star_name: str,
                           ensemble_vmag: Optional[float] = None) -> int:
    """Write a simple photometry CSV with differential and optional apparent magnitudes.

    Columns: JD, V_C (differential V-minus-ensemble), V_app (if ensemble_vmag known), err.
    All rows included (including high-error ones); NaN magnitudes skipped.
    Returns number of data rows written.
    """
    rows_written = 0
    with open(dat) as f_in, open(out, "w", newline="", encoding="utf-8") as f_out:
        f_out.write(f"# Star: {star_name}\n")
        if ensemble_vmag is not None:
            f_out.write(f"# Ensemble V (APASS comp stars median): {ensemble_vmag:.3f}\n")
            f_out.write("# V_app = V_C + ensemble_V  (approximate apparent V magnitude)\n")
        w = csv.writer(f_out)
        if ensemble_vmag is not None:
            w.writerow(["JD", "V_C", "V_app", "err"])
        else:
            w.writerow(["JD", "V_C", "err"])

        for line in f_in:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                jd_i = next(i for i, p in enumerate(parts) if float(p) > 2_400_000)
                jd   = float(parts[jd_i])
                vc   = float(parts[jd_i + 1])
                err  = float(parts[jd_i + 2]) if len(parts) > jd_i + 2 else 0.0
                if math.isnan(vc):
                    continue
                if ensemble_vmag is not None:
                    w.writerow([f"{jd:.6f}", f"{vc:.4f}",
                                f"{vc + ensemble_vmag:.4f}", f"{err:.4f}"])
                else:
                    w.writerow([f"{jd:.6f}", f"{vc:.4f}", f"{err:.4f}"])
                rows_written += 1
            except (ValueError, StopIteration, IndexError):
                continue
    return rows_written


def generate_light_curve_plot(dat: Path, out: Path, star_name: str,
                               ensemble_vmag: Optional[float] = None) -> bool:
    """Generate a two-panel matplotlib light curve (differential + apparent if known).

    Returns True if the plot was saved successfully.
    """
    if not HAS_MATPLOTLIB:
        return False

    jds, vcs, errs = [], [], []
    with open(dat) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                jd_i = next(i for i, p in enumerate(parts) if float(p) > 2_400_000)
                jd   = float(parts[jd_i])
                vc   = float(parts[jd_i + 1])
                err  = float(parts[jd_i + 2]) if len(parts) > jd_i + 2 else 0.0
                if math.isnan(vc) or math.isnan(err):
                    continue
                jds.append(jd); vcs.append(vc); errs.append(err)
            except (ValueError, StopIteration, IndexError):
                continue

    if not jds:
        return False

    jd0 = int(min(jds))
    xs  = [j - jd0 for j in jds]

    has_apparent = ensemble_vmag is not None
    n_panels     = 2 if has_apparent else 1
    fig, axes    = plt.subplots(n_panels, 1,
                                figsize=(11, 4 * n_panels),
                                squeeze=False,
                                sharex=True)
    fig.suptitle(star_name, fontsize=13, fontweight="bold")

    kw = dict(fmt="o", ms=3, elinewidth=0.8, capsize=2)

    ax1 = axes[0][0]
    ax1.errorbar(xs, vcs, yerr=errs, color="#5294e2", ecolor="#aaaaaa", **kw)
    ax1.set_ylabel("V − C  (differential)")
    ax1.invert_yaxis()
    ax1.grid(alpha=0.25)
    ax1.set_title("Differential photometry  (V − C)")

    if has_apparent:
        vapps = [v + ensemble_vmag for v in vcs]
        ax2 = axes[1][0]
        ax2.errorbar(xs, vapps, yerr=errs, color="#7ab648", ecolor="#aaaaaa", **kw)
        ax2.set_ylabel("V  (apparent)")
        ax2.invert_yaxis()
        ax2.grid(alpha=0.25)
        ax2.set_title(f"Apparent magnitude  "
                      f"(ensemble comp V = {ensemble_vmag:.2f})")
        ax2.set_xlabel(f"JD − {jd0}")
    else:
        ax1.set_xlabel(f"JD − {jd0}")

    plt.tight_layout()
    try:
        plt.savefig(out, dpi=150, bbox_inches="tight")
    except Exception:
        plt.close(fig)
        return False
    plt.close(fig)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────────────────────────────────────

def run_pipeline(config: dict,
                 on_log: Callable[[str], None],
                 on_progress: Callable[[int, str], None],
                 on_done: Callable[[bool, str], None]) -> None:
    """Run the full photometry pipeline.

    config keys:
      session    (Path)             session root folder
      star       (dict)             {name, ra, dec, mag, …}
      comp_stars (list[dict])       fallback comparison stars from APASS
      nstars     (int)              max comparison stars for findcompstars
      start_step (int)              1, 3, 4, or 5
      dark_dir   (Optional[Path])
      flat_dir   (Optional[Path])
      bias_dir   (Optional[Path])
      filt_code  (str)              AAVSO filter code, e.g. "CV"
      filt_note  (str)              extra note, e.g. "LP_filter_Seestar_S30Pro"
      runner     (SirilRunner)
      phot_aperture  (float)        setphot forced-aperture radius px (def 10)
      phot_inner     (float)        sky-annulus inner radius px (def 20)
      phot_outer     (float)        sky-annulus outer radius px (def 30)
      phot_dyn_ratio (float)        dynamic aperture = 0.5*FWHM*ratio (def 4.0)
      phot_min_val   (float)        min valid pixel value (def -1000; registered
                                    frames are background-subtracted so sky sits
                                    slightly below 0 on the 16-bit ADU scale)
      phot_max_val   (float)        max valid pixel value / saturation (def 60000)
      error_floor    (bool)         add a systematic error floor to MERR so it
                                    reflects the real scatter, not just Siril's
                                    formal photon noise (def True)
      error_floor_mag (Optional[float]) fixed σ_syst (mag) to add in quadrature
                                    instead of measuring it from the light curve

    Callbacks on_log/on_progress/on_done are called from the pipeline thread.
    """
    session    = config["session"]
    star       = config["star"]
    comp_stars = list(config.get("comp_stars", []))
    nstars     = config.get("nstars", 10)
    start_step = config.get("start_step", 1)
    filt_code  = config.get("filt_code", "CV")
    filt_note  = config.get("filt_note", "")
    runner     = config["runner"]

    # ── Photometry (setphot) parameters — user-tunable ────────────────────────
    # See Step 5 for how these map onto Siril's aperture-photometry model.
    # Defaults reproduce the historical hard-coded values.
    phot = {
        "aperture":  float(config.get("phot_aperture",  10.0)),
        "inner":     float(config.get("phot_inner",     20.0)),
        "outer":     float(config.get("phot_outer",     30.0)),
        "dyn_ratio": float(config.get("phot_dyn_ratio",  4.0)),
        "min_val":   float(config.get("phot_min_val", -1000.0)),
        "max_val":   float(config.get("phot_max_val", 60000.0)),
    }

    proc    = session / "process"
    masters = proc / "masters"
    proc.mkdir(exist_ok=True)

    # ── Resolve the lights directory (tolerant to layout) ─────────────────────
    #   legacy Seestar dump : session/lights/*.fits
    #   Argos               : session/Lights/<filter>/*.fits
    lights = resolve_frame_dir(session, ["lights", "light"]) or (session / "lights")

    # ── Calibration frames: honour an explicit dir, else auto-detect ──────────
    def _resolve_calib(cfg_dir, names):
        if cfg_dir is not None:
            p = Path(cfg_dir)
            if _dir_has_fits(p):
                return p
            # user pointed at a container (e.g. Darks/) — descend to the frames
            subs = [s for s in sorted(p.iterdir()) if _dir_has_fits(s)] if p.is_dir() else []
            if subs:
                return max(subs, key=lambda s: len(list_fits(s)))
            return p  # leave as-is; the calibration step will report the failure
        return resolve_frame_dir(session, names)

    dark_dir = _resolve_calib(config.get("dark_dir"), ["darks", "dark"])
    flat_dir = _resolve_calib(config.get("flat_dir"), ["flats", "flat"])
    bias_dir = _resolve_calib(config.get("bias_dir"), ["biases", "bias"])

    # ── Per-star output subfolder (results/<safe_name>/) ──────────────────────
    safe        = star["name"].replace(" ", "_").replace("/", "-")
    results_dir = session / "results" / safe
    results_dir.mkdir(parents=True, exist_ok=True)

    # ── Log file — accumulate entries; write atomically when pipeline ends ─────
    _log: list[str] = [
        f"# Seestar Variable Star Pipeline v{VERSION}",
        f"# Target : {star['name']}  RA={star['ra']:.5f}  Dec={star['dec']:+.5f}",
        f"# Session: {session}",
        f"# Started: {datetime.datetime.now().isoformat()}",
        "",
    ]
    _orig_on_log = on_log

    def on_log(msg: str) -> None:  # shadows the callback parameter
        _log.append(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}")
        _orig_on_log(msg)

    _orig_on_done = on_done

    def on_done(ok: bool, msg: str) -> None:  # shadows the callback parameter
        _log.append(f"\n[{datetime.datetime.now().strftime('%H:%M:%S')}] "
                    f"{'SUCCESS' if ok else 'FAILED'}: {msg}")
        try:
            (results_dir / "pipeline.log").write_text(
                "\n".join(_log) + "\n", encoding="utf-8")
        except Exception:
            pass
        _orig_on_done(ok, msg)

    # ── Report the resolved input directories ─────────────────────────────────
    if lights.is_dir():
        on_log(f"Lights: {lights.relative_to(session) if session in lights.parents else lights}"
               f"  ({len(list_fits(lights))} frames)")
    else:
        on_log(f"WARNING: lights directory not found under {session}")
    for label, d in [("Darks", dark_dir), ("Flats", flat_dir), ("Bias", bias_dir)]:
        if d is not None and _dir_has_fits(d):
            on_log(f"{label}: {d}  ({len(list_fits(d))} frames)")

    # Warn on dark/light exposure mismatch — dark subtraction needs matching t_exp.
    if dark_dir is not None and _dir_has_fits(dark_dir):
        t_light = fits_exptime(lights) if lights.is_dir() else None
        t_dark  = fits_exptime(dark_dir)
        if t_light and t_dark and abs(t_light - t_dark) > 0.5:
            on_log(f"WARNING: dark exposure {t_dark:g}s ≠ light exposure {t_light:g}s — "
                   "dark subtraction may be inaccurate.")

    focal = SEESTAR["focal"]
    pixsz = SEESTAR["pixsz"]
    gain  = get_gain_eadu(lights)
    on_log(f"Gain: {gain} e-/ADU")

    # ── Jump to step 3/4/5 — detect existing registered sequence ─────────────
    if start_step == 3:
        for candidate in ("pp_light_", "light_"):
            if (proc / f"{candidate}.seq").exists():
                seq = candidate
                break
        else:
            on_done(False,
                    f"No registered sequence (light_.seq) found in {proc.name}/.\n"
                    "Run the full pipeline first (steps 1–5).")
            return
        registered = f"r_{seq}"

    if start_step >= 4:
        for candidate in ("r_pp_light_", "r_light_"):
            if (proc / f"{candidate}.seq").exists():
                registered = candidate
                break
        else:
            on_done(False,
                    f"No registered sequence found in {proc.name}/.\n"
                    "Run the full pipeline first (steps 1–5).")
            return
        on_log(f"Resuming from step {start_step} — sequence: {registered}")

    # ── Steps 0–3 (only when starting from step 1) ───────────────────────────
    if start_step == 1:
        calib_dirs = {
            k: d for k, d in [("dark", dark_dir), ("flat", flat_dir), ("bias", bias_dir)]
            if d is not None and _dir_has_fits(d)
        }
        has_calib = bool(calib_dirs)

        bad_out = lights / f'"{proc}"'
        if bad_out.exists():
            shutil.rmtree(bad_out, ignore_errors=True)
            on_log("Cleaned up misplaced output directory from a previous run.")

        # ── Step 1: link raw lights → light_ sequence ─────────────────────────
        # Must run before calibration: `calibrate` operates on the light_ sequence.
        on_progress(15, "Linking frames…")
        on_log("─── Step 1: link lights → sequence ───")
        if not runner.run_script(proc, [
            f'cd "{lights}"',
            f'link light -out={proc}',
        ], "_s1_link.ssf"):
            on_done(False, "Step 1 failed (link lights)")
            return

        # ── Step 0: build calibration masters & calibrate light_ → pp_light_ ──
        if has_calib:
            masters.mkdir(exist_ok=True)
            on_progress(22, "Building calibration masters…")
            on_log("─── Step 0: calibration frames ───")
            calib_cmds = []
            for kind, src in calib_dirs.items():
                norm = "-nonorm" if kind in ("dark", "bias") else "-norm=mul"
                calib_cmds += [
                    f'cd "{src}"',
                    f'link {kind} -out={masters}',
                    f'cd "{masters}"',
                    f'stack {kind}_ rej 3 3 {norm} -out=master_{kind}',
                ]
                on_log(f"  → stacking {kind}s from {src}/")
            cal_flags = " ".join(
                f"-{k}=masters/master_{k}"
                for k in ("bias", "dark", "flat") if k in calib_dirs
            )
            # Cosmetic correction detects hot/cold pixels from the master dark;
            # it requires -dark, so only enable it when a dark is present.
            cc = " -cc=dark" if "dark" in calib_dirs else ""
            calib_cmds += [f'cd "{proc}"', f"calibrate light_ {cal_flags}{cc}"]
            if not runner.run_script(proc, calib_cmds, "_s0_calibrate.ssf"):
                on_done(False, "Step 0 failed (calibration)")
                return
            seq = "pp_light_"
        else:
            on_log("No calibration frames — proceeding with raw lights.")
            seq = "light_"

        # ── Step 2: register the (calibrated) sequence ────────────────────────
        on_progress(30, "Computing registration…")
        on_log("─── Step 2: register -2pass ───")
        if not runner.run_script(proc, [
            f'cd "{proc}"',
            f"register {seq} -2pass",
        ], "_s2_register.ssf"):
            on_done(False, "Step 2 failed (register)")
            return

        on_progress(45, "Aligning frames…")
        on_log("─── Step 3: seqapplyreg ───")
        if not runner.run_script(proc, [
            f'cd "{proc}"',
            f"seqapplyreg {seq} -framing=min -filter-round=2.5k",
        ], "_s3_applyreg.ssf"):
            on_done(False, "Step 3 failed (seqapplyreg)")
            return

        registered = f"r_{seq}"

    if start_step == 3:
        on_progress(45, "Aligning frames…")
        on_log("─── Step 3: seqapplyreg ───")
        if not runner.run_script(proc, [
            f'cd "{proc}"',
            f"seqapplyreg {seq} -framing=min -filter-round=2.5k",
        ], "_s3_applyreg.ssf"):
            on_done(False, "Step 3 failed (seqapplyreg)")
            return

    # ── Step 4: plate solve ───────────────────────────────────────────────────
    if start_step <= 4:
        on_progress(65, "Plate solving registered frames…")
        on_log("─── Step 4: seqplatesolve ───")
        disto = proc / "ps_distortion"
        disto_arg = "-disto=ps_distortion" if disto.is_dir() else ""
        if not runner.run_script(proc, [
            f'cd "{proc}"',
            f"seqplatesolve {registered} -nocache -force "
            f"-focal={focal} -pixelsize={pixsz} -radius=2.5 {disto_arg}".strip(),
        ], "_s4_platesolve.ssf"):
            on_log("WARNING: plate solve had errors — continuing")

    # ── Step 5: sequence integrity + photometry ───────────────────────────────
    on_progress(82, "Checking sequence integrity…")

    seq_file    = proc / f"{registered}.seq"
    stem        = registered.rstrip("_")
    seq_fixlen  = 4
    seq_ref_idx = None
    img_entries: list[tuple[int, int]] = []
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

    if img_entries:
        valid_nums = {n for n, _ in img_entries}
        for f in sorted(proc.glob(f"{stem}_*.fit")):
            try:
                num = int(f.stem[len(stem)+1:])
                if num not in valid_nums:
                    on_log(f"Removing stale frame {f.name}")
                    f.unlink()
            except (ValueError, OSError):
                pass

    ref_img_num = None
    if seq_ref_idx is not None and seq_ref_idx < len(img_entries):
        ref_img_num = img_entries[seq_ref_idx][0]

    on_progress(85, "Aperture photometry…")
    on_log("─── Step 5: setphot + light_curve ───")

    comp_csv = proc / "comp_stars.csv"
    lc_cmd   = None

    # ── Resolve ref frame WCS once (reused by frame check + FALLBACK A) ──────
    ref_hdr: dict = {}
    naxis1 = naxis2 = 0
    if ref_img_num is not None:
        _rf = proc / f"{stem}_{ref_img_num:0{seq_fixlen}d}.fit"
        if _rf.exists():
            ref_hdr = read_fits_header(_rf)
            naxis1  = int(ref_hdr.get("NAXIS1", 0))
            naxis2  = int(ref_hdr.get("NAXIS2", 0))

    def _to_disp(fx: float, fy: float) -> tuple[int, int]:
        """WCS pixel (1-indexed, y-up FITS) → Siril display pixel (0-indexed, y-down)."""
        return round(fx - 0.5), round(naxis2 - fy + 0.5)

    # ── Early frame check (before any findcompstars / APASS query) ────────────
    # For alt-az mounts, field rotation limits the reliable area to the inscribed
    # circle of the frame (radius = min(NAXIS1, NAXIS2) / 2).  Stars outside this
    # circle may fall into the black rotation corners of some registered frames.
    # With -framing=min, black corners are already eliminated, but the inscribed
    # circle check is still the correct conservative bound for the VSX-filtered
    # candidates.
    tdx: Optional[int] = None
    tdy: Optional[int] = None
    SAFE_MARGIN = 50   # px buffer inside inscribed circle (absorbs drift + aperture)
    if naxis1 and naxis2:
        tx, ty = sky_to_pixel(star["ra"], star["dec"], ref_hdr)
        if tx is not None:
            tdx, tdy = _to_disp(tx, ty)
            # Inscribed circle check in display pixel space
            cx, cy = naxis1 / 2.0, naxis2 / 2.0
            safe_r  = min(naxis1, naxis2) / 2.0 - SAFE_MARGIN
            dist_sq = (tdx - cx) ** 2 + (tdy - cy) ** 2
            if dist_sq > safe_r ** 2:
                dist_px = dist_sq ** 0.5
                on_done(False,
                        f"'{star['name']}' is outside the alt-az safe zone.\n"
                        f"Distance from field centre: {dist_px:.0f} px, "
                        f"safe radius: {safe_r:.0f} px "
                        f"(inscribed circle of {naxis1}×{naxis2} − {SAFE_MARGIN} px margin).\n"
                        "Field rotation on an alt-az mount means this star lands in the "
                        "black corners of many registered frames.\n"
                        "Re-run the VSX query and choose a star closer to the field centre.")
                return
            on_log(f"Target at display pixel ({tdx}, {tdy}) — inside alt-az safe zone ✓")
        else:
            on_log("Frame WCS not available — skipping bounds check")

    ensemble_vmag: Optional[float] = None   # comp-star median V; enables apparent mag

    # ── PRIMARY: findcompstars → -ninastars ───────────────────────────────────
    # VSX names sometimes have a "V* " prefix that GCVS/SIMBAD doesn't use.
    ref_fits_name = (
        f"{stem}_{ref_img_num:0{seq_fixlen}d}.fit"
        if ref_img_num is not None else None
    )
    if ref_fits_name and (proc / ref_fits_name).exists():
        if comp_csv.exists():
            comp_csv.unlink()
        star_arg = star["name"].removeprefix("V* ").replace('"', '\\"')
        on_log(f"[PRIMARY] findcompstars '{star_arg}' (APASS, dvmag=3, emag=0.05)")
        runner.run_script(proc, [
            f'cd "{proc}"',
            f"load {ref_fits_name}",
            f'findcompstars "{star_arg}" -narrow -dvmag=3 -emag=0.05 -catalog=apass -out=comp_stars.csv',
        ], "_s5a_findcomp.ssf")
        if comp_csv.exists() and comp_csv.stat().st_size > 50:
            kept = truncate_comp_csv(comp_csv, max(3, min(50, nstars)))
            ensemble_vmag = _parse_comp_csv_vmag(comp_csv)
            if ensemble_vmag is not None:
                on_log(f"[PRIMARY] ensemble comp V = {ensemble_vmag:.2f} "
                       f"(median of {kept} APASS stars from comp_stars.csv)")
            lc_cmd = f"light_curve {registered} 0 -ninastars=comp_stars.csv"
            on_log(f"[PRIMARY] OK — {kept} comp stars via findcompstars")
        else:
            on_log(f"[PRIMARY] FAILED — '{star_arg}' not in Siril catalog → trying fallback A")

    # ── FALLBACK A: APASS via VizieR + pixel coords ───────────────────────────
    # dvmag=3.5: VSX reports max brightness, star may be 1-2 mag fainter at minimum.
    if lc_cmd is None and ref_img_num is not None:
        if not comp_stars:
            on_log(f"[FALLBACK A] querying APASS (radius=1.5°, dvmag=3.5) …")
            comp_stars = query_apass(star["ra"], star["dec"], star["mag"],
                                     radius_deg=1.5, n=nstars, dvmag=3.5)
            if comp_stars:
                on_log(f"[FALLBACK A] {len(comp_stars)} APASS comp stars found")
            else:
                on_log(f"[FALLBACK A] no APASS stars (target mag={star['mag']:.1f}) → trying fallback B")

    if lc_cmd is None and comp_stars and tdx is not None and naxis1:
        margin = 35
        ref_pix: list[tuple[int, int]] = []
        _ens_comps: list[dict] = []   # in-frame comps tracked for ensemble Vmag
        for cs in comp_stars:
            rx, ry = sky_to_pixel(cs["ra"], cs["dec"], ref_hdr)
            if rx is not None:
                rdx, rdy = _to_disp(rx, ry)
                if (margin < rdx < naxis1 - margin and
                        margin < rdy < naxis2 - margin):
                    ref_pix.append((rdx, rdy))
                    _ens_comps.append(cs)
        if not ref_pix and ref_hdr.get("CRVAL1"):
            # Target near the field edge — retry APASS around the frame centre so
            # at least some reference stars land inside the image.
            on_log("[FALLBACK A] no in-frame comp stars near target "
                   "— retrying APASS around field centre (radius=0.5°) …")
            centre_comps = query_apass(
                float(ref_hdr["CRVAL1"]), float(ref_hdr["CRVAL2"]),
                star["mag"], radius_deg=0.5, n=nstars, dvmag=3.5)
            for cs in centre_comps:
                rx, ry = sky_to_pixel(cs["ra"], cs["dec"], ref_hdr)
                if rx is not None:
                    rdx, rdy = _to_disp(rx, ry)
                    if (margin < rdx < naxis1 - margin and
                            margin < rdy < naxis2 - margin):
                        ref_pix.append((rdx, rdy))
                        _ens_comps.append(cs)
            if ref_pix:
                on_log(f"[FALLBACK A] {len(ref_pix)} comp stars from field centre")
            else:
                on_log("[FALLBACK A] no comp stars in frame (tried target + field centre)")
        if ref_pix:
            _vmags = [cs["vmag"] for cs in _ens_comps if cs.get("vmag", 0) > 0]
            if _vmags:
                ensemble_vmag = statistics.median(_vmags)
                on_log(f"[FALLBACK A] ensemble comp V = {ensemble_vmag:.2f} "
                       f"(median of {len(_vmags)} APASS stars)")
            lc_cmd = f"light_curve {registered} 0 -at={tdx},{tdy}"
            for rdx, rdy in ref_pix:
                lc_cmd += f" -refat={rdx},{rdy}"
            on_log(f"[FALLBACK A] OK — target ({tdx},{tdy}), "
                   f"{len(ref_pix)} comp stars in frame")
        else:
            on_log("[FALLBACK A] no comp stars in frame → trying fallback B")

    # ── FALLBACK B: sky coordinates (-wcs/-refwcs) ────────────────────────────
    # Only used when we have no pixel coords (no WCS on ref frame).
    # Always filter comp stars to in-frame only — Siril aborts on the first
    # -refwcs coordinate that falls outside the image.
    if lc_cmd is None:
        if not comp_stars:
            on_done(False,
                    f"No comparison stars found for '{star['name']}'.\n"
                    "Check: internet connection for APASS query, target magnitude "
                    "filter, or select a different star.")
            return
        inframe_comps: list[dict] = []
        if tdx is not None and naxis1:
            margin = 35
            for cs in comp_stars:
                rx, ry = sky_to_pixel(cs["ra"], cs["dec"], ref_hdr)
                if rx is not None:
                    rdx, rdy = _to_disp(rx, ry)
                    if (margin < rdx < naxis1 - margin and
                            margin < rdy < naxis2 - margin):
                        inframe_comps.append(cs)
        else:
            inframe_comps = comp_stars  # no WCS — pass all and let Siril decide
        if not inframe_comps:
            on_done(False,
                    f"No in-frame comparison stars found for '{star['name']}'.\n"
                    "The target may be at the edge of the field with no suitable "
                    "APASS reference stars in the image. Try a different star.")
            return
        _vmags = [cs["vmag"] for cs in inframe_comps if cs.get("vmag", 0) > 0]
        if _vmags:
            ensemble_vmag = statistics.median(_vmags)
        on_log("[FALLBACK B] using sky coordinates (-wcs/-refwcs)")
        lc_cmd = (f"light_curve {registered} 0 "
                  f"-wcs={star['ra']:.6f},{star['dec']:.6f}")
        for cs in inframe_comps:
            lc_cmd += f" -refwcs={cs['ra']:.6f},{cs['dec']:.6f}"

    # Delete any stale light_curve.dat from a previous star's run so we never
    # silently pick up the wrong data if Siril fails to write a new one.
    stale_dat = proc / "light_curve.dat"
    if stale_dat.exists():
        stale_dat.unlink()
        on_log("Removed stale light_curve.dat from previous run")

    # Strip BAYERPAT before light_curve to work around a Siril 1.4.3 bug:
    # copyfits(CP_FORMAT) clears date_obs in the CFA copy used for PSF fitting,
    # so seq->imgparam[i].date_obs is never populated → JD axis shows frame indices.
    # Without BAYERPAT, Siril skips the CFA copy and reads date_obs correctly.
    bayer_stripped: dict[Path, tuple[int, str]] = {}
    for fit_path in sorted(proc.glob(f"{stem}_*.fit")):
        result = _fits_strip_keyword(fit_path, "BAYERPAT")
        if result:
            bayer_stripped[fit_path] = result
    if bayer_stripped:
        on_log(f"BAYERPAT stripped from {len(bayer_stripped)} frames (date_obs fix)")

    # Dynamic aperture = 0.5 * FWHM * dyn_ratio (Siril photometry.c). A frame's
    # PSF fit fails with "inner radii too small (N required)" when that dynamic
    # aperture reaches `inner`, i.e. FWHM >= 2*inner/dyn_ratio px. Raise inner/
    # outer (or lower dyn_ratio) to keep marginal frames; raise max_val when the
    # comparison stars saturate ("reference stars ... out of valid pixel range").
    setphot_cmd = (
        f"setphot -aperture={phot['aperture']:g} -inner={phot['inner']:g} "
        f"-outer={phot['outer']:g} -dyn_ratio={phot['dyn_ratio']:g} "
        f"-min_val={phot['min_val']:g} -max_val={phot['max_val']:g} -gain={gain}"
    )
    on_log(f"Photometry: aperture={phot['aperture']:g} inner={phot['inner']:g} "
           f"outer={phot['outer']:g} dyn_ratio={phot['dyn_ratio']:g} "
           f"valid=[{phot['min_val']:g},{phot['max_val']:g}]")
    ok = runner.run_script(proc, [
        f'cd "{proc}"',
        setphot_cmd,
        lc_cmd,
    ], "_s5_phot.ssf")

    for fit_path, (offset, val) in bayer_stripped.items():
        _fits_restore_keyword(fit_path, "BAYERPAT", offset, val)
    if bayer_stripped:
        on_log("BAYERPAT restored")

    # Siril light_curve exits non-zero when ANY frame fails PSF fitting (e.g. target
    # lands in a black corner on some frames).  It still writes light_curve.dat for
    # the frames that succeeded.  Treat a non-empty dat file as partial success.
    lc_dat = proc / "light_curve.dat"
    if not ok:
        if lc_dat.exists() and lc_dat.stat().st_size > 50:
            on_log("Step 5 partially succeeded — some frames failed PSF, "
                   "processing available data")
        else:
            on_done(False, "Step 5 failed (light_curve — no output produced)")
            return

    # ── Collect results ───────────────────────────────────────────────────────
    # results_dir / safe defined at top of function; subfolder already created.
    out_dat  = results_dir / "light_curve.dat"
    out_csv  = results_dir / "aavso.csv"
    out_phot = results_dir / "photometry.csv"
    out_png  = results_dir / "light_curve.png"

    if lc_dat.exists():
        try:
            first_fit = next(iter(sorted(proc.glob(f"{stem}_*.fit"))))
            exptime_s = float(read_fits_header(first_fit).get("EXPTIME", 20.0))
        except StopIteration:
            exptime_s = 20.0
        if inject_jd_into_dat(lc_dat, seq_file, stem, seq_fixlen, proc, exptime_s):
            on_log("JD normalized to absolute JD ✓")
        shutil.copy2(lc_dat, out_dat)

        # Realistic uncertainties: Siril's formal error underestimates the true
        # scatter (flat/scintillation/calibration systematics absent from the photon
        # budget). Lift MERR to sqrt(formal^2 + sigma_syst^2) from the measured scatter.
        # All downstream exports read out_dat, so correcting it once is enough.
        if config.get("error_floor", True):
            ef = apply_error_floor(out_dat, config.get("error_floor_mag"))
            if ef and ef["sigma_syst"] > 0:
                on_log(f"Error floor: measured scatter {ef['sigma_real']:.4f} mag vs "
                       f"formal median {ef['median_err']:.4f} → added σ_syst="
                       f"{ef['sigma_syst']:.4f} mag in quadrature to MERR")
            elif ef:
                on_log("Error floor: formal errors already consistent with scatter — "
                       "MERR unchanged")
            else:
                on_log("Error floor: too few points to measure scatter — MERR unchanged")

        # AAVSO extended format (differential V-C, ENSEMBLE comp)
        # Airmass per point comes from the registered frame headers (AIRMASS keyword)
        airmass_map = build_airmass_map(proc, stem, seq_fixlen)
        if airmass_map:
            on_log(f"Airmass: read from {len(airmass_map)} frame headers")
        else:
            on_log("Airmass: no AIRMASS keyword in frame headers — AMASS left as 'na'")
        rows = export_aavso(out_dat, out_csv, star["name"],
                            filt_code=filt_code, filt_note=filt_note,
                            airmass_map=airmass_map)
        n_valid = len(rows)
        n_total = sum(1 for ln in out_dat.read_text().splitlines()
                      if ln and not ln.startswith("#"))
        on_log(f"AAVSO CSV: {n_valid}/{n_total} pts exported "
               f"(NaN and MERR>0.5 excluded)")

        # Simple photometry CSV with optional apparent magnitude
        n_phot = export_photometry_csv(out_dat, out_phot, star["name"], ensemble_vmag)
        if ensemble_vmag is not None:
            on_log(f"Photometry CSV: {n_phot} pts, apparent V = V_C + {ensemble_vmag:.2f}")
        else:
            on_log(f"Photometry CSV: {n_phot} pts (apparent mag unavailable — PRIMARY path)")

        # FWHM per frame from registration .seq (plate scale ~1.035"/px for Seestar)
        out_fwhm = results_dir / "fwhm.csv"
        n_fwhm = export_fwhm_csv(seq_file, proc, stem, seq_fixlen, out_fwhm,
                                  plate_scale_arcsec=1.035)
        if n_fwhm:
            on_log(f"FWHM CSV: {n_fwhm} frames written → {out_fwhm.name}")
        else:
            on_log("FWHM CSV: skipped (no registration FWHM in .seq)")

        # Matplotlib plot (two panels if apparent mag known; fallback to Siril's plot)
        plotted = generate_light_curve_plot(out_dat, out_png, star["name"], ensemble_vmag)
        if not plotted:
            png = proc / "light_curve.png"
            if png.exists():
                shutil.copy2(png, out_png)

        on_done(True, str(results_dir))
    else:
        on_done(False,
                "light_curve.dat not found.\n"
                "The target may be outside the field, or too faint.")
