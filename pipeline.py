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
import subprocess
from pathlib import Path
from typing import Callable, Optional

VERSION = "0.1.1"

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
        x = math.degrees(-math.cos(dec_r) * math.sin(dra) / denom)
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
            if header_skipped <= 2:
                continue
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


def get_gain_eadu(lights_dir: Path) -> float:
    """Read real gain in e-/ADU from the first FITS file in lights/.
    GAIN=200 on a Seestar is ISO-equivalent, not e-/ADU.
    Returns 1.0 as a safe fallback.
    """
    for fit in sorted(lights_dir.glob("*.fit"))[:1]:
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

def export_aavso(dat: Path, out: Path, name: str,
                 filt_code: str = "CV", filt_note: str = "",
                 merr_max: float = 0.5) -> list:
    """Export light curve to AAVSO Extended format.

    Filters NaN magnitudes and MERR > merr_max.
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
            w.writerow([name, f"{jd:.6f}", f"{mag:.4f}", f"{err:.4f}",
                        filt_code, "NO", "DIFF",
                        "ENSEMBLE", "na", "na", "na", "na", "1", "na", notes])
    return rows


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

    Callbacks on_log/on_progress/on_done are called from the pipeline thread.
    """
    session    = config["session"]
    star       = config["star"]
    comp_stars = list(config.get("comp_stars", []))
    nstars     = config.get("nstars", 10)
    start_step = config.get("start_step", 1)
    dark_dir   = config.get("dark_dir")
    flat_dir   = config.get("flat_dir")
    bias_dir   = config.get("bias_dir")
    filt_code  = config.get("filt_code", "CV")
    filt_note  = config.get("filt_note", "")
    runner     = config["runner"]

    proc    = session / "process"
    lights  = session / "lights"
    masters = proc / "masters"
    proc.mkdir(exist_ok=True)

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
        has_calib = any([dark_dir, flat_dir, bias_dir])

        bad_out = lights / f'"{proc}"'
        if bad_out.exists():
            shutil.rmtree(bad_out, ignore_errors=True)
            on_log("Cleaned up misplaced output directory from a previous run.")

        if has_calib:
            masters.mkdir(exist_ok=True)
            on_progress(5, "Building calibration masters…")
            on_log("─── Step 0: calibration frames ───")
            calib_cmds = []
            for kind, src in [("dark", dark_dir), ("flat", flat_dir), ("bias", bias_dir)]:
                if src is None:
                    continue
                norm = "-nonorm" if kind in ("dark", "bias") else "-norm=mul"
                calib_cmds += [
                    f'cd "{src}"',
                    f'link {kind} -out={masters}',
                    f'cd "{masters}"',
                    f'stack {kind}_ rej 3 3 {norm} -out=master_{kind}',
                ]
                on_log(f"  → stacking {kind}s from {src.name}/")
            cal_flags = " ".join(
                f"-{k}=masters/master_{k}"
                for k, d in [("bias", bias_dir), ("dark", dark_dir), ("flat", flat_dir)]
                if d is not None
            )
            calib_cmds += [f'cd "{proc}"', f"calibrate light_ {cal_flags} -cc=banding"]
            if not runner.run_script(proc, calib_cmds, "_s0_calibrate.ssf"):
                on_done(False, "Step 0 failed (calibration)")
                return
            seq = "pp_light_"
        else:
            on_log("No calibration frames — proceeding with raw lights.")
            seq = "light_"

        on_progress(15, "Linking frames & computing registration…")
        on_log("─── Step 1: link lights → sequence ───")
        on_log("─── Step 2: register -2pass ───")
        if not runner.run_script(proc, [
            f'cd "{lights}"',
            f'link light -out={proc}',
            f'cd "{proc}"',
            f"register {seq} -2pass",
        ], "_s12_link_register.ssf"):
            on_done(False, "Step 1/2 failed (link / register)")
            return

        on_progress(45, "Aligning frames…")
        on_log("─── Step 3: seqapplyreg ───")
        if not runner.run_script(proc, [
            f'cd "{proc}"',
            f"seqapplyreg {seq} -framing=cog -filter-round=2.5k",
        ], "_s3_applyreg.ssf"):
            on_done(False, "Step 3 failed (seqapplyreg)")
            return

        registered = f"r_{seq}"

    if start_step == 3:
        on_progress(45, "Aligning frames…")
        on_log("─── Step 3: seqapplyreg ───")
        if not runner.run_script(proc, [
            f'cd "{proc}"',
            f"seqapplyreg {seq} -framing=cog -filter-round=2.5k",
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

    # ── PRIMARY: findcompstars → -ninastars ───────────────────────────────────
    # Requires a plate-solved ref frame loaded into gfit (has_wcs check in Siril).
    # VSX names sometimes have a "V* " prefix that GCVS doesn't use — strip it.
    ref_fits_name = (
        f"{stem}_{ref_img_num:0{seq_fixlen}d}.fit"
        if ref_img_num is not None else None
    )
    if ref_fits_name and (proc / ref_fits_name).exists():
        if comp_csv.exists():
            comp_csv.unlink()
        # Normalize name: strip "V* " prefix that VSX adds but GCVS/Siril doesn't use
        star_arg = star["name"].removeprefix("V* ").replace('"', '\\"')
        on_log(f"[PRIMARY] findcompstars '{star_arg}' (APASS, dvmag=3, emag=0.05)")
        runner.run_script(proc, [
            f'cd "{proc}"',
            f"load {ref_fits_name}",
            f'findcompstars "{star_arg}" -narrow -dvmag=3 -emag=0.05 -catalog=apass -out=comp_stars.csv',
        ], "_s5a_findcomp.ssf")
        if comp_csv.exists() and comp_csv.stat().st_size > 50:
            kept = truncate_comp_csv(comp_csv, max(3, min(50, nstars)))
            lc_cmd = f"light_curve {registered} 0 -ninastars=comp_stars.csv"
            on_log(f"[PRIMARY] OK — {kept} comp stars via findcompstars")
        else:
            on_log(f"[PRIMARY] FAILED — '{star_arg}' not found in Siril catalog "
                   f"(GCVS/SIMBAD name mismatch?) → trying fallback A")

    # ── FALLBACK A: APASS via VizieR + pixel coords from WCS ─────────────────
    # Uses dvmag=3.5 because VSX reports max brightness, not mean — the star may
    # be 1-2 mag fainter than listed when observed near minimum.
    if lc_cmd is None and ref_img_num is not None:
        if not comp_stars:
            on_log(f"[FALLBACK A] querying APASS (radius=1.5°, dvmag=3.5) …")
            comp_stars = query_apass(star["ra"], star["dec"], star["mag"],
                                     radius_deg=1.5, n=nstars, dvmag=3.5)
            if comp_stars:
                on_log(f"[FALLBACK A] {len(comp_stars)} APASS comp stars found")
            else:
                on_log(f"[FALLBACK A] no APASS stars found "
                       f"(target mag={star['mag']:.1f}, radius=1.5°) → trying fallback B")

    if lc_cmd is None and ref_img_num is not None and comp_stars:
        ref_fits = proc / f"{stem}_{ref_img_num:0{seq_fixlen}d}.fit"
        ref_hdr  = read_fits_header(ref_fits)
        naxis1   = int(ref_hdr.get("NAXIS1", 0))
        naxis2   = int(ref_hdr.get("NAXIS2", 0))
        tx, ty   = sky_to_pixel(star["ra"], star["dec"], ref_hdr)

        if tx is None or naxis1 == 0 or naxis2 == 0:
            on_log("[FALLBACK A] WCS not available on ref frame → trying fallback B")
        else:
            def _to_disp(fx, fy, n2=naxis2):
                return round(fx - 0.5), round(n2 - fy + 0.5)

            tdx, tdy = _to_disp(tx, ty)

            # Margin = outer ring + small buffer so PSF fit doesn't touch the edge
            margin = 35
            if not (margin < tdx < naxis1 - margin and
                    margin < tdy < naxis2 - margin):
                on_done(False,
                        f"'{star['name']}' is outside the image frame "
                        f"(display pos: {tdx},{tdy}; image: {naxis1}×{naxis2}).\n"
                        "Select a star closer to the center of the field.")
                return

            ref_pix = []
            for cs in comp_stars:
                rx, ry = sky_to_pixel(cs["ra"], cs["dec"], ref_hdr)
                if rx is not None:
                    rdx, rdy = _to_disp(rx, ry)
                    if (margin < rdx < naxis1 - margin and
                            margin < rdy < naxis2 - margin):
                        ref_pix.append((rdx, rdy))

            if ref_pix:
                lc_cmd = f"light_curve {registered} 0 -at={tdx},{tdy}"
                for rdx, rdy in ref_pix:
                    lc_cmd += f" -refat={rdx},{rdy}"
                on_log(f"[FALLBACK A] OK — target ({tdx},{tdy}), "
                       f"{len(ref_pix)}/{len(comp_stars)} comp stars in frame")
            else:
                on_log("[FALLBACK A] no comp stars project inside frame "
                       "→ trying fallback B")

    # ── FALLBACK B: sky coordinates (-wcs/-refwcs) ────────────────────────────
    # Last resort — less reliable, but avoids a silent failure.
    if lc_cmd is None:
        if not comp_stars:
            on_done(False,
                    f"No comparison stars found for '{star['name']}'.\n"
                    "Check: internet connection for APASS query, target magnitude "
                    "filter, or select a different star.")
            return
        on_log("[FALLBACK B] using sky coordinates (-wcs/-refwcs)")
        lc_cmd = (f"light_curve {registered} 0 "
                  f"-wcs={star['ra']:.6f},{star['dec']:.6f}")
        for cs in comp_stars:
            lc_cmd += f" -refwcs={cs['ra']:.6f},{cs['dec']:.6f}"

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

    ok = runner.run_script(proc, [
        f'cd "{proc}"',
        f"setphot -aperture=10 -inner=20 -outer=30 -dyn_ratio=4.0 -gain={gain}",
        lc_cmd,
    ], "_s5_phot.ssf")

    for fit_path, (offset, val) in bayer_stripped.items():
        _fits_restore_keyword(fit_path, "BAYERPAT", offset, val)
    if bayer_stripped:
        on_log("BAYERPAT restored")

    if not ok:
        on_done(False, "Step 5 failed (light_curve)")
        return

    # ── Collect results ───────────────────────────────────────────────────────
    lc_dat = proc / "light_curve.dat"
    results_dir = session / "results"
    results_dir.mkdir(exist_ok=True)
    safe = star["name"].replace(" ", "_").replace("/", "-")
    out_dat = results_dir / f"{safe}_light_curve.dat"
    out_csv = results_dir / f"{safe}_aavso.csv"

    if lc_dat.exists():
        try:
            first_fit = next(iter(sorted(proc.glob(f"{stem}_*.fit"))))
            exptime_s = float(read_fits_header(first_fit).get("EXPTIME", 20.0))
        except StopIteration:
            exptime_s = 20.0
        if inject_jd_into_dat(lc_dat, seq_file, stem, seq_fixlen, proc, exptime_s):
            on_log("JD normalized to absolute JD ✓")
        shutil.copy2(lc_dat, out_dat)

        out_png = results_dir / f"{safe}_light_curve.png"
        png = proc / "light_curve.png"
        if png.exists():
            shutil.copy2(png, out_png)

        rows = export_aavso(out_dat, out_csv, star["name"],
                            filt_code=filt_code, filt_note=filt_note)
        n_valid = len(rows)
        n_total = sum(1 for l in out_dat.read_text().splitlines()
                      if l and not l.startswith("#"))
        on_log(f"AAVSO CSV: {n_valid}/{n_total} pts exported "
               f"(NaN and MERR>0.5 excluded)")
        on_done(True, str(out_dat))
    else:
        on_done(False,
                "light_curve.dat not found.\n"
                "The target may be outside the field, or too faint.")
