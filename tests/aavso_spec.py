"""Checker for the AAVSO Extended File Format, written from the specification.

https://www.aavso.org/aavso-extended-file-format (version 1.2, last updated
2021-12-02). Kept independent of pipeline.py on purpose: the exporter is tested
against the spec, not against itself.
"""

from pathlib import Path

REQUIRED_PARAMS = ("TYPE", "OBSCODE", "SOFTWARE", "DELIM", "DATE")
FIELDS = ("NAME", "DATE", "MAG", "MERR", "FILT", "TRANS", "MTYPE", "CNAME",
          "CMAG", "KNAME", "KMAG", "AMASS", "GROUP", "CHART", "NOTES")

# FILTER list of the spec, minus "O" which WebObs currently rejects.
FILTERS = {
    "U", "B", "V", "R", "I", "J", "H", "K", "TG", "TB", "TR", "CV", "CR",
    "SZ", "SU", "SG", "SR", "SI", "STU", "STV", "STB", "STY", "STHBW", "STHBN",
    "MA", "MB", "MI", "ZS", "Y", "HA", "HAC",
}


def _is_number(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def _check_optional_number(problems, where, field, value, limit):
    if value.lower() != "na" and not _is_number(value):
        problems.append(f"{where}: {field} {value!r} is neither a number nor 'na'")
    if len(value) > limit:
        problems.append(f"{where}: {field} longer than {limit} characters")


def aavso_violations(path) -> list[str]:
    """Every rule of the spec that the file at *path* breaks (empty if none)."""
    path = Path(path)
    problems: list[str] = []
    if path.suffix.lower() not in (".txt", ".csv", ".tsv"):
        problems.append("extension must be .txt, .csv or .tsv")

    lines = path.read_text(encoding="utf-8").split("\n")
    if lines and lines[-1] == "":
        lines.pop()

    params: dict[str, str] = {}
    data: list[tuple[int, str]] = []
    for n, line in enumerate(lines, start=1):
        if not line.strip():
            problems.append(f"line {n}: blank line not commented out with #")
        elif line.startswith("#"):
            key, sep, value = line[1:].partition("=")
            if sep and key.strip().upper() in (*REQUIRED_PARAMS, "OBSTYPE"):
                if data and key.strip().upper() not in ("OBSCODE", "DATE"):
                    problems.append(f"line {n}: #{key} after the first observation")
                params[key.strip().upper()] = value.strip()
        else:
            data.append((n, line))

    for p in REQUIRED_PARAMS:
        if p not in params:
            problems.append(f"missing #{p}= parameter")
    if params.get("TYPE", "").upper() != "EXTENDED":
        problems.append("#TYPE must be EXTENDED")
    if "OBSCODE" in params and not params["OBSCODE"]:
        problems.append("#OBSCODE is empty")
    if len(params.get("SOFTWARE", "")) > 255:
        problems.append("#SOFTWARE longer than 255 characters")
    delim = {"tab": "\t", "comma": ","}.get(params.get("DELIM", "").lower(),
                                          params.get("DELIM", ""))
    if len(delim) != 1 or (delim != "\t" and not 32 <= ord(delim) <= 126) \
            or delim in "|# ":
        problems.append(f"#DELIM {params.get('DELIM')!r} is not an allowed delimiter")
        return problems
    date_fmt = params.get("DATE", "").upper()
    if date_fmt not in ("JD", "HJD", "EXCEL"):
        problems.append(f"#DATE {params.get('DATE')!r} is not JD, HJD or EXCEL")
    if params.get("OBSTYPE", "CCD").upper() not in ("CCD", "DSLR", "PEP"):
        problems.append(f"#OBSTYPE {params['OBSTYPE']!r} is not CCD, DSLR or PEP")
    if not data:
        problems.append("no observation lines")

    for n, line in data:
        where = f"line {n}"
        values = line.split(delim)
        if len(values) != len(FIELDS):
            problems.append(f"{where}: {len(values)} fields instead of {len(FIELDS)}")
            continue
        for field, value in zip(FIELDS, values):
            if value != value.strip():
                problems.append(f"{where}: {field} has a leading or trailing space")
            if not value:
                problems.append(f"{where}: {field} is empty (use 'na')")
        f = {k: v.strip() for k, v in zip(FIELDS, values)}

        if len(f["NAME"]) > 30:
            problems.append(f"{where}: NAME longer than 30 characters")
        if date_fmt in ("JD", "HJD") and not _is_number(f["DATE"]):
            problems.append(f"{where}: DATE {f['DATE']!r} is not a Julian date")
        if len(f["DATE"]) > 16:
            problems.append(f"{where}: DATE longer than 16 characters")
        mag = f["MAG"].removeprefix("<")
        if not _is_number(mag):
            problems.append(f"{where}: MAG {f['MAG']!r} is not a magnitude")
        if "." not in mag:
            problems.append(f"{where}: MAG {f['MAG']!r} has no decimal point")
        if len(f["MAG"]) > 8:
            problems.append(f"{where}: MAG longer than 8 characters")
        _check_optional_number(problems, where, "MERR", f["MERR"], 6)
        if f["FILT"].upper() not in FILTERS:
            problems.append(f"{where}: FILT {f['FILT']!r} is not an AAVSO filter")
        if f["TRANS"].upper() not in ("YES", "NO"):
            problems.append(f"{where}: TRANS {f['TRANS']!r} is not YES or NO")
        mtype = f["MTYPE"].upper()
        if mtype not in ("STD", "DIF", "ABS"):
            problems.append(f"{where}: MTYPE {f['MTYPE']!r} is not STD or DIF")
        if len(f["CNAME"]) > 20:
            problems.append(f"{where}: CNAME longer than 20 characters")
        _check_optional_number(problems, where, "CMAG", f["CMAG"], 8)
        if f["CNAME"].upper() == "ENSEMBLE" and f["CMAG"].lower() != "na":
            problems.append(f"{where}: ensemble photometry needs CMAG=na")
        if mtype == "DIF" and (f["CNAME"].lower() == "na"
                               or f["CNAME"].upper() == "ENSEMBLE"
                               or not _is_number(f["CMAG"])):
            problems.append(f"{where}: MTYPE=DIF needs a comparison star and its CMAG")
        if len(f["KNAME"]) > 20:
            problems.append(f"{where}: KNAME longer than 20 characters")
        _check_optional_number(problems, where, "KMAG", f["KMAG"], 8)
        _check_optional_number(problems, where, "AMASS", f["AMASS"], 7)
        if f["GROUP"].lower() != "na" and not f["GROUP"].isdigit():
            problems.append(f"{where}: GROUP {f['GROUP']!r} is neither an integer nor 'na'")
        if len(f["GROUP"]) > 5:
            problems.append(f"{where}: GROUP longer than 5 characters")
        if len(f["CHART"]) > 20:
            problems.append(f"{where}: CHART longer than 20 characters")
    return problems
