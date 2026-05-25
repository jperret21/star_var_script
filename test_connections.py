#!/usr/bin/env python3
"""
Connection & API diagnostics for seestar_varstar_siril.py
Run this first if you see catalog or network errors.

Usage:
    python3 test_connections.py
    # or from Siril: Script > Run Script > test_connections.py
"""

import sys
import socket
import json

# ── Try requests ──────────────────────────────────────────────────────────────
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# M81 field coordinates used as test target
TEST_RA  = 149.28
TEST_DEC = 69.09

PASS = "✓  PASS"
FAIL = "✗  FAIL"
WARN = "⚠  WARN"

results = []

def check(name, ok, detail=""):
    tag = PASS if ok else FAIL
    line = f"  {tag}  {name}"
    if detail:
        line += f"\n         {detail}"
    print(line)
    results.append(ok)
    return ok


# ─────────────────────────────────────────────────────────────────────────────
print("=" * 60)
print("  Seestar Variable Star — Connection Diagnostics")
print("=" * 60)

# ── 1. Python & requests ──────────────────────────────────────────────────────
print("\n[1] Python environment")
check("Python version", True, f"{sys.version.split()[0]}")
check("requests module", HAS_REQUESTS,
      "pip install requests" if not HAS_REQUESTS else "")

# ── 2. Basic internet ─────────────────────────────────────────────────────────
print("\n[2] Internet connectivity")
try:
    socket.setdefaulttimeout(5)
    socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect(("8.8.8.8", 53))
    check("Internet (DNS reachable)", True)
except Exception as e:
    check("Internet (DNS reachable)", False, str(e))

# ── 3. VizieR VSX (primary variable star catalog) ─────────────────────────────
print("\n[3] VizieR / VSX  (primary variable star source)")
if HAS_REQUESTS:
    try:
        url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
        params = {
            "-source": "B/vsx/vsx",
            "-c": f"{TEST_RA} {TEST_DEC}",
            "-c.r": "30", "-c.u": "arcmin",
            "-out": "Name,Type,RAJ2000,DEJ2000",
            "-out.max": "3",
        }
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        rows = [l for l in r.text.splitlines()
                if l.strip() and not l.startswith("#") and not l.startswith("-")]
        # rows[0] = column headers, rows[1] = units, rows[2+] = data
        n_stars = max(0, len(rows) - 2)
        check("VizieR VSX query", True,
              f"{n_stars} stars returned (M81 field, r=30 arcmin)")
        if n_stars > 0 and len(rows) > 2:
            print(f"         Sample: {rows[2][:70].strip()}")
    except Exception as e:
        check("VizieR VSX query", False, str(e))
else:
    check("VizieR VSX query", False, "requests not available")

# ── 4. VizieR APASS (comparison stars) ───────────────────────────────────────
print("\n[4] VizieR / APASS  (comparison stars)")
if HAS_REQUESTS:
    try:
        url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
        params = {
            "-source": "II/336/apass9",
            "-c": f"{TEST_RA} {TEST_DEC}",
            "-c.r": "20", "-c.u": "arcmin",
            "-out": "RAJ2000,DEJ2000,Vmag",
            "-out.max": "5",
        }
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        rows = [l for l in r.text.splitlines()
                if l.strip() and not l.startswith("#") and not l.startswith("-")]
        n_stars = max(0, len(rows) - 2)
        check("VizieR APASS query", True,
              f"{n_stars} comparison stars (M81 field, r=20 arcmin)")
    except Exception as e:
        check("VizieR APASS query", False, str(e))
else:
    check("VizieR APASS query", False, "requests not available")

# ── 5. AAVSO VSX API (secondary, legacy) ─────────────────────────────────────
print("\n[5] AAVSO VSX API  (secondary source — may be unreliable)")
if HAS_REQUESTS:
    try:
        url = "https://vsx.aavso.org/index.php"
        params = {"view": "api.list", "ra": TEST_RA, "dec": TEST_DEC,
                  "radius": "30", "format": "json"}
        r = requests.get(url, params=params, timeout=15, allow_redirects=True)
        r.raise_for_status()
        data = r.json()
        objs = data.get("VSXObjects", {}).get("VSXObject", [])
        if not isinstance(objs, list):
            objs = [objs] if objs else []
        check("AAVSO VSX API", True, f"{len(objs)} objects returned")
    except json.JSONDecodeError as e:
        check("AAVSO VSX API", False,
              f"JSON parse error — API may be down. VizieR VSX will be used instead.")
    except Exception as e:
        check("AAVSO VSX API", False, f"{e} — VizieR VSX will be used instead.")
else:
    check("AAVSO VSX API", False, "requests not available")

# ── 6. sirilpy ───────────────────────────────────────────────────────────────
print("\n[6] Siril integration")
try:
    sys.path.insert(0, "/Applications/Siril.app/Contents/Resources/share/siril/python_module")
    import sirilpy
    from sirilpy import SirilInterface
    iface = SirilInterface()
    iface.connect()
    check("sirilpy connect", True, "Siril is running and connected")
    try:
        wd = iface.get_wd()
        check("Siril working directory", True, str(wd))
    except Exception:
        check("Siril working directory", True, "(get_wd not available — use Browse)")
except ImportError:
    check("sirilpy module", False, "Run this script from Siril: Script > Run Script")
except Exception as e:
    check("sirilpy connect", False,
          f"{e}\nOpen Siril first, then run this script from Script > Run Script")

# ── 7. siril-cli ──────────────────────────────────────────────────────────────
import shutil
cli = None
for c in ["siril-cli", "/Applications/Siril.app/Contents/MacOS/siril-cli"]:
    if shutil.which(c):
        cli = c
        break
check("siril-cli executable", cli is not None,
      cli if cli else "Install Siril 1.4+ from https://siril.org")

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
passed = sum(results)
total  = len(results)
color_ok = passed == total
print(f"  Result: {passed}/{total} checks passed")
if not color_ok:
    print("  Fix the FAIL items above before running the main script.")
else:
    print("  All good — you can run seestar_varstar_siril.py")
print("=" * 60)
