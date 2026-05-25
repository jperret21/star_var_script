"""Tests for AAVSO export and comparison star CSV truncation."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import csv
import shutil
import tempfile
import unittest
from pathlib import Path

from pipeline import export_aavso, truncate_comp_csv, VERSION


def _write_dat(path: Path, rows: list[str]) -> Path:
    lines = ["#JD_UT", "#Star"] + rows
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_comp_csv(path: Path, n_comp: int) -> Path:
    lines = ["Target,ES UMa,148.888,69.065,12.5"]
    for i in range(1, n_comp + 1):
        lines.append(f"Comp1,star{i},{148+i*0.01:.3f},69.0,{12+i*0.1:.1f},0.02,12.0")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class TestExportAavso(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _read_csv_rows(self, path: Path):
        rows = []
        with open(path, newline="") as f:
            for line in f:
                if line.startswith("#"):
                    continue
            f.seek(0)
            reader = csv.DictReader(
                (l for l in f if not l.startswith("#"))
            )
            for row in reader:
                rows.append(row)
        return rows

    def _read_headers(self, path: Path):
        return [l.strip() for l in path.read_text().splitlines()
                if l.startswith("#")]

    def test_valid_rows_written(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.0230",
            "2460325.439000 12.3670 0.0310",
        ])
        out = self.tmp / "out.csv"
        rows = export_aavso(dat, out, "ES UMa")
        self.assertEqual(len(rows), 2)
        self.assertTrue(out.exists())

    def test_nan_magnitude_excluded(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.0230",
            "2460325.439000 nan 0.0",
        ])
        out = self.tmp / "out.csv"
        rows = export_aavso(dat, out, "ES UMa")
        self.assertEqual(len(rows), 1)

    def test_high_merr_excluded(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.0230",
            "2460325.439000 12.3670 0.6000",  # MERR > 0.5
        ])
        out = self.tmp / "out.csv"
        rows = export_aavso(dat, out, "ES UMa")
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(rows[0][2], 0.0230, places=4)

    def test_merr_at_threshold_excluded(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.5001",
        ])
        rows = export_aavso(dat, self.tmp / "out.csv", "ES UMa")
        self.assertEqual(len(rows), 0)

    def test_merr_at_threshold_included(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.4999",
        ])
        rows = export_aavso(dat, self.tmp / "out.csv", "ES UMa")
        self.assertEqual(len(rows), 1)

    def test_empty_dat_returns_empty(self):
        dat = _write_dat(self.tmp / "lc.dat", [])
        rows = export_aavso(dat, self.tmp / "out.csv", "ES UMa")
        self.assertEqual(rows, [])
        self.assertFalse((self.tmp / "out.csv").exists())

    def test_csv_headers_present(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        export_aavso(dat, out, "ES UMa")
        headers = self._read_headers(out)
        self.assertIn("#TYPE=EXTENDED", headers)
        self.assertIn("#DELIM=,", headers)
        self.assertIn("#DATE=JD", headers)
        self.assertTrue(any("SOFTWARE" in h for h in headers))
        self.assertTrue(any(VERSION in h for h in headers))

    def test_filter_code_in_csv(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        export_aavso(dat, out, "ES UMa", filt_code="V")
        headers = self._read_headers(out)
        self.assertIn("#FILTER=V", headers)

    def test_star_name_in_data_rows(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        export_aavso(dat, out, "ES UMa")
        rows = self._read_csv_rows(out)
        self.assertEqual(rows[0]["NAME"], "ES UMa")

    def test_filt_note_in_notes_column(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        export_aavso(dat, out, "ES UMa", filt_note="LP_filter")
        rows = self._read_csv_rows(out)
        self.assertIn("LP_filter", rows[0]["NOTES"])

    def test_custom_merr_max(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.8000",
        ])
        rows = export_aavso(dat, self.tmp / "out.csv", "ES UMa", merr_max=1.0)
        self.assertEqual(len(rows), 1)


class TestTruncateCompCsv(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_truncates_to_n(self):
        f = _write_comp_csv(self.tmp / "comp.csv", n_comp=10)
        kept = truncate_comp_csv(f, 5)
        self.assertEqual(kept, 5)
        lines = [l for l in f.read_text().splitlines()
                 if l.startswith("Comp1,")]
        self.assertEqual(len(lines), 5)

    def test_keeps_target_line(self):
        f = _write_comp_csv(self.tmp / "comp.csv", n_comp=5)
        truncate_comp_csv(f, 3)
        lines = [l for l in f.read_text().splitlines()
                 if l.startswith("Target,")]
        self.assertEqual(len(lines), 1)

    def test_fewer_than_n_keeps_all(self):
        f = _write_comp_csv(self.tmp / "comp.csv", n_comp=3)
        kept = truncate_comp_csv(f, 10)
        self.assertEqual(kept, 3)

    def test_n_zero_removes_all_comp(self):
        f = _write_comp_csv(self.tmp / "comp.csv", n_comp=5)
        kept = truncate_comp_csv(f, 0)
        self.assertEqual(kept, 0)
        lines = [l for l in f.read_text().splitlines()
                 if l.startswith("Comp1,")]
        self.assertEqual(len(lines), 0)


if __name__ == "__main__":
    unittest.main()
