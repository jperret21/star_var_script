"""Shared helpers for pipeline tests — no external dependencies."""

from pathlib import Path


def make_fits(path: Path, keywords: dict) -> Path:
    """Write a minimal valid FITS file with the given header keywords.

    keywords: {name: value}  — strings stored as FITS string cards,
                               numbers stored as FITS numeric cards.
    """
    def card(kw: str, val) -> bytes:
        kw8 = kw.upper().ljust(8)[:8]
        if isinstance(val, str):
            field = f"'{val:<8}'"
        elif isinstance(val, bool):
            field = "T" if val else "F"
            field = field.rjust(20)
        else:
            field = str(val).rjust(20)
        raw = f"{kw8}= {field}"
        return raw.encode("ascii").ljust(80)[:80]

    cards = [card("SIMPLE", True)]
    for k, v in keywords.items():
        cards.append(card(k, v))
    end = b"END" + b" " * 77
    cards.append(end)

    raw = b"".join(cards)
    rem = len(raw) % 2880
    if rem:
        raw += b" " * (2880 - rem)

    path.write_bytes(raw)
    return path


def make_seq(path: Path, selected: list[int], all_count: int) -> Path:
    """Write a minimal Siril .seq file.

    selected: list of 1-based image numbers that are selected (flag=1).
    all_count: total number of I-lines to write (1..all_count).
    """
    lines = [f"I {n} {1 if n in selected else 0}" for n in range(1, all_count + 1)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def make_dat(path: Path, header_line: str, rows: list[str]) -> Path:
    """Write a minimal Siril light_curve.dat file."""
    lines = [header_line, "#Frame magnitude error"] + rows
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
