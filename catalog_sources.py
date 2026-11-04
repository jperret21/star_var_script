"""
catalog_sources.py — variable-star and comparison-star catalogue queries.

Each catalogue has a primary source (VizieR) and a fallback (the AAVSO APIs):

  variable stars     VizieR B/vsx/vsx       →  AAVSO VSX  (vsx.aavso.org api.list)
  comparison stars   VizieR II/336/apass9   →  AAVSO VSP  (app.aavso.org/vsp/api/chart/)

A source is tried only when the previous one times out or fails; an empty but
valid answer is kept. Network only, no Siril dependency: it can be tried from a
terminal (python3 catalog_sources.py RA DEC).

AAVSO API references:
  VSX  https://archive.aavso.org/apis-aavso-resources  (api.list: ra, dec,
       radius in degrees, tomag; format=json)
  VSP  https://archive.aavso.org/api-vsp  (ra, dec, fov in arcmin, maglimit)
"""

from __future__ import annotations

import math
import re
import time
from typing import Callable, Optional, Sequence

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

VIZIER_URL    = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
# www.aavso.org/vsx is behind a Cloudflare browser check; vsx.aavso.org serves
# the same API without it.
AAVSO_VSX_URL = "https://vsx.aavso.org/index.php"
AAVSO_VSP_URL = "https://app.aavso.org/vsp/api/chart/"

# Order in which the sources are tried ("vizier", "aavso").
SOURCE_ORDER: tuple[str, ...] = ("vizier", "aavso")
TIMEOUT = 20          # s, per request
# A source that timed out or failed is skipped for this long, so that one
# pipeline run does not wait for the same dead server at every query.
RETRY_AFTER = 600     # s

VSX_MAX_ROWS = 500    # as VizieR -out.max: the nearest stars are kept
# Above a 180′ field VSP accepts no maglimit fainter than 12 (HTTP 400).
VSP_MAX_FOV = 180.0   # arcmin

APASS_LABEL = "APASS DR9"
VSP_LABEL   = "AAVSO VSP"

SOURCE_NAMES = {
    ("vsx", "vizier"):   "VizieR VSX",
    ("vsx", "aavso"):    "AAVSO VSX",
    ("comps", "vizier"): "VizieR APASS",
    ("comps", "aavso"):  "AAVSO VSP",
}

Log = Optional[Callable[[str], None]]


class CatalogError(Exception):
    """A source did not answer, or answered something unreadable."""


_failed_at: dict[str, float] = {}


def source_down(source: str) -> bool:
    """True while *source* is skipped after a recent failure."""
    t = _failed_at.get(source)
    return t is not None and time.monotonic() - t < RETRY_AFTER


def reset_failures() -> None:
    _failed_at.clear()


def _get(url: str, params: dict) -> "requests.Response":
    if not HAS_REQUESTS:
        raise CatalogError("the requests module is not available")
    try:
        r = requests.get(url, params=params, timeout=TIMEOUT)
        r.raise_for_status()
    except requests.Timeout:
        raise CatalogError(f"no answer within {TIMEOUT} s") from None
    except requests.HTTPError as e:
        raise CatalogError(f"HTTP {e.response.status_code}") from None
    except requests.RequestException as e:
        raise CatalogError(type(e).__name__) from None
    return r


def _json(r: "requests.Response"):
    try:
        return r.json()
    except ValueError:
        raise CatalogError("the answer is not JSON") from None


def _sep_deg(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    """Angular separation in degrees."""
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    c = (math.sin(d1) * math.sin(d2)
         + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _sexa_to_deg(text: str, hours: bool) -> float:
    """'02:46:49.70' (hours) or '+41:10:11.6' (degrees) → degrees."""
    s = text.strip()
    sign = -1.0 if s.startswith("-") else 1.0
    parts = [float(p) for p in s.lstrip("+-").split(":")]
    value = parts[0] + (parts[1] if len(parts) > 1 else 0) / 60 \
        + (parts[2] if len(parts) > 2 else 0) / 3600
    return sign * value * (15.0 if hours else 1.0)


# ─────────────────────────────────────────────────────────────────────────────
# VizieR
# ─────────────────────────────────────────────────────────────────────────────

def vizier_tsv(catalog: str, ra: float, dec: float,
               radius_arcmin: float, columns: list[str],
               filters: dict | None = None, max_rows: int = 500) -> list[list[str]]:
    """Data rows of a VizieR cone search. Raises CatalogError."""
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
    r = _get(VIZIER_URL, params)
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
    if not in_data and "<html" in r.text[:500].lower():
        raise CatalogError("the answer is not a VizieR table")
    return rows


def vizier_vsx(ra: float, dec: float, radius_deg: float) -> list[dict]:
    """Variable stars from VizieR B/vsx/vsx."""
    rows = vizier_tsv(
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


def vizier_apass(ra: float, dec: float, target_mag: float,
                 radius_deg: float, dvmag: float) -> list[dict]:
    """APASS DR9 stars from VizieR II/336/apass9, unsorted.

    The magnitude range is filtered by VizieR: rows come sorted by distance
    and capped at max_rows, so filtering here instead would only see the few
    hundred nearest (mostly faint) stars — for a bright target, none in range.
    """
    rows = vizier_tsv(
        "II/336/apass9", ra, dec, radius_deg * 60,
        ["RAJ2000", "DEJ2000", "Vmag", "e_Vmag"],
        filters={"e_Vmag": "<0.05",
                 "Vmag": f"{target_mag - dvmag:.2f}..{target_mag + dvmag:.2f}"},
        max_rows=200,
    )
    comps = []
    for parts in rows:
        if len(parts) < 3:
            continue
        try:
            comps.append({"ra": float(parts[0].strip()),
                          "dec": float(parts[1].strip()),
                          "vmag": float(parts[2].strip()),
                          "source": APASS_LABEL, "chart": APASS_LABEL})
        except (ValueError, IndexError):
            continue
    return comps


# ─────────────────────────────────────────────────────────────────────────────
# AAVSO
# ─────────────────────────────────────────────────────────────────────────────

_LEADING_NUMBER = re.compile(r"^\s*[<>]?\s*(-?\d+(?:\.\d+)?)")


def _vsx_mag(text: Optional[str]) -> float:
    """'10.8 V', '<15.2 G', '13.1: r' → leading number, 99.0 if none."""
    m = _LEADING_NUMBER.match(text or "")
    return float(m.group(1)) if m else 99.0


def aavso_vsx(ra: float, dec: float, radius_deg: float) -> list[dict]:
    """Variable stars from the AAVSO VSX api.list, the VSX_MAX_ROWS nearest."""
    data = _json(_get(AAVSO_VSX_URL, {
        "view": "api.list", "ra": f"{ra:.6f}", "dec": f"{dec:+.6f}",
        "radius": f"{radius_deg:.6f}", "format": "json",
    }))
    try:
        objs = (data.get("VSXObjects") or {}).get("VSXObject") or []
    except AttributeError:
        raise CatalogError("unexpected VSX answer") from None
    if isinstance(objs, dict):       # a single star is not wrapped in a list
        objs = [objs]
    stars = []
    for o in objs:
        try:
            s_ra, s_dec = float(o["RA2000"]), float(o["Declination2000"])
        except (KeyError, TypeError, ValueError):
            continue
        period = o.get("Period")
        try:
            period_s = f"{float(period):.4f} d" if period else "—"
        except ValueError:
            period_s = "—"
        stars.append({
            "name":     (o.get("Name") or "?").strip(),
            "ra":       s_ra,
            "dec":      s_dec,
            "mag":      _vsx_mag(o.get("MaxMag")),
            "var_type": (o.get("VariabilityType") or "").strip() or "?",
            "period":   period_s,
            "dist":     0.0,
            "auid":     o.get("AUID") or "",
        })
    stars.sort(key=lambda s: _sep_deg(ra, dec, s["ra"], s["dec"]))
    return stars[:VSX_MAX_ROWS]


def aavso_vsp(ra: float, dec: float, target_mag: float,
              radius_deg: float, dvmag: float) -> list[dict]:
    """Sequence stars of an AAVSO VSP chart centred on (ra, dec), unsorted.

    The chart is a square of side 2·radius, at most VSP_MAX_FOV; only the stars
    within radius_deg are kept, as for a cone search. Stars without a V
    magnitude are left out; e_vmag is None where VSP gives no error.
    """
    fov = min(2 * radius_deg * 60, VSP_MAX_FOV)
    maglimit = math.ceil(round((target_mag + dvmag) * 10, 6)) / 10
    data = _json(_get(AAVSO_VSP_URL, {
        "format": "json", "ra": f"{ra:.6f}", "dec": f"{dec:+.6f}",
        "fov": f"{fov:g}", "maglimit": f"{maglimit:g}",
    }))
    if not isinstance(data, dict) or "photometry" not in data:
        raise CatalogError("unexpected VSP answer")
    chart = data.get("chartid") or None
    comps = []
    for p in data.get("photometry") or []:
        try:
            c_ra  = _sexa_to_deg(str(p["ra"]), hours=True)
            c_dec = _sexa_to_deg(str(p["dec"]), hours=False)
        except (KeyError, ValueError):
            continue
        v = next((b for b in p.get("bands") or [] if b.get("band") == "V"), None)
        if v is None or v.get("mag") is None:
            continue
        if _sep_deg(ra, dec, c_ra, c_dec) > radius_deg:
            continue
        comps.append({"ra": c_ra, "dec": c_dec, "vmag": float(v["mag"]),
                      "e_vmag": v.get("error"),
                      "auid": p.get("auid") or "", "label": str(p.get("label") or ""),
                      "source": VSP_LABEL, "chart": chart})
    return comps


# ─────────────────────────────────────────────────────────────────────────────
# Fallback
# ─────────────────────────────────────────────────────────────────────────────

def _first_answer(kind: str, calls: dict[str, Callable[[], list[dict]]],
                  order: Optional[Sequence[str]], on_log: Log
                  ) -> tuple[list[dict], Optional[str]]:
    """Rows of the first source in *order* that answers, and its name."""
    log = on_log or (lambda _msg: None)
    for source in (order or SOURCE_ORDER):
        if source not in calls:
            log(f"Catalogue: unknown source '{source}' ignored")
            continue
        name = SOURCE_NAMES[(kind, source)]
        if source_down(source):
            log(f"{name}: skipped (failed less than {RETRY_AFTER // 60} min ago)")
            continue
        try:
            rows = calls[source]()
        except CatalogError as e:
            _failed_at[source] = time.monotonic()
            log(f"{name}: failed ({e})")
            continue
        _failed_at.pop(source, None)
        log(f"{name}: {len(rows)} stars")
        return rows, name
    log("Catalogue: no source answered")
    return [], None


def fetch_vsx(ra: float, dec: float, radius_deg: float, *,
              order: Optional[Sequence[str]] = None,
              on_log: Log = None) -> tuple[list[dict], Optional[str]]:
    """Variable stars within radius_deg, and the name of the source used."""
    return _first_answer("vsx", {
        "vizier": lambda: vizier_vsx(ra, dec, radius_deg),
        "aavso":  lambda: aavso_vsx(ra, dec, radius_deg),
    }, order, on_log)


def fetch_comps(ra: float, dec: float, target_mag: float, *,
                radius_deg: float = 1.0, n: int = 6, dvmag: float = 2.0,
                order: Optional[Sequence[str]] = None,
                on_log: Log = None) -> tuple[list[dict], Optional[str]]:
    """Comparison stars closest in V to target_mag, and the source used.

    Kept: |V − target_mag| ≤ dvmag, V error < 0.05 where known, not within 5″
    of (ra, dec). Each star has ra, dec, vmag, plus source (APASS_LABEL or
    VSP_LABEL), chart (the VSP chart ID, or APASS_LABEL) and, from VSP, auid.
    """
    rows, name = _first_answer("comps", {
        "vizier": lambda: vizier_apass(ra, dec, target_mag, radius_deg, dvmag),
        "aavso":  lambda: aavso_vsp(ra, dec, target_mag, radius_deg, dvmag),
    }, order, on_log)
    comps = []
    for c in rows:
        if abs(c["vmag"] - target_mag) > dvmag:
            continue
        err = c.get("e_vmag")
        if err is not None and err >= 0.05:
            continue
        sep = ((c["ra"] - ra) * 3600) ** 2 + ((c["dec"] - dec) * 3600) ** 2
        if sep < 25:  # < 5 arcsec from target
            continue
        comps.append(c)
    comps.sort(key=lambda c: abs(c["vmag"] - target_mag))
    return comps[:n], name


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Try the catalogue sources from a terminal.")
    ap.add_argument("ra", type=float, help="degrees")
    ap.add_argument("dec", type=float, help="degrees")
    ap.add_argument("--radius", type=float, default=1.5, help="degrees (default 1.5)")
    ap.add_argument("--mag", type=float, default=11.0, help="target V (default 11)")
    ap.add_argument("--dvmag", type=float, default=3.5)
    ap.add_argument("--order", default=",".join(SOURCE_ORDER),
                    help="comma-separated sources (default vizier,aavso)")
    args = ap.parse_args()
    order = [s.strip() for s in args.order.split(",") if s.strip()]
    stars, src = fetch_vsx(args.ra, args.dec, args.radius, order=order, on_log=print)
    for s in stars[:10]:
        print(f"  {s['name']:<30} {s['ra']:10.5f} {s['dec']:+10.5f}  "
              f"max {s['mag']:5.2f}  {s['var_type']:<8} {s['period']}")
    comps, src = fetch_comps(args.ra, args.dec, args.mag, radius_deg=args.radius,
                             n=10, dvmag=args.dvmag, order=order, on_log=print)
    for c in comps:
        print(f"  {c.get('auid', ''):<12} {c['ra']:10.5f} {c['dec']:+10.5f}  "
              f"V {c['vmag']:6.3f}  {c['source']} {c['chart'] or ''}")
