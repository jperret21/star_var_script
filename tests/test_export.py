"""Tests for AAVSO export and comparison star CSV truncation."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import csv
import math
import shutil
import statistics
import tempfile
import unittest
from pathlib import Path

from pipeline import (
    FILTER_OPTIONS,
    VERSION,
    _parse_comp_csv,
    ensemble_zero_point,
    export_aavso,
    siril_loaded_comps,
    siril_refs_used,
    truncate_comp_csv,
)
from tests.aavso_spec import FILTERS, aavso_violations


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


def _export(dat: Path, out: Path, name: str = "ES UMa", **kw) -> list:
    """export_aavso with a valid observer code and ensemble magnitude."""
    kw.setdefault("obscode", "TST01")
    kw.setdefault("ensemble_vmag", 12.0)
    kw.setdefault("n_comps", 5)
    return export_aavso(dat, out, name, **kw)


# Examples 1–3 of the spec page (trailing spaces left by its HTML removed).
_SPEC_EXAMPLES = [
    """#TYPE=EXTENDED
#OBSCODE=TST01
#SOFTWARE=IRAF 12.4
#DELIM=,
#DATE=JD
#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
SS CYG,2450702.1234,11.135,0.003,V,NO,STD,105,na,na,na,na,na,X16382L,This is a test
""",
    """#TYPE=EXTENDED
#OBSCODE=TST01
#SOFTWARE=GCX 2.0
#DELIM=,
#DATE=JD
#OBSTYPE=CCD
#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
SS CYG,2450702.1234,11.235,0.003,B,NO,STD,105,10.593,110,11.090,1.561,1,X16382L,na
SS CYG,2450702.1254,11.135,0.003,V,NO,STD,105,10.594,110,10.994,1.563,1,X16382L,na
""",
    """#TYPE=EXTENDED
#OBSCODE=TST01
#SOFTWARE=GCX 2.0
#DELIM=,
#DATE=JD
#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
SS CYG,2450702.1234,11.235,0.003,B,NO,STD,ENSEMBLE,na,105,10.593,1.561,1,X16382L,na
SS CYG,2450702.1254,11.135,0.003,V,NO,STD,ENSEMBLE,na,105,10.492,1.563,1,X16382L,na
""",
]

# What v0.1.5 wrote — rejected by WebObs with
# "invalid literal for int() with base 10: np.str_('MTYPE')".
_OLD_EXPORT = """#TYPE=EXTENDED
#OBSCODE=XXXX
#SOFTWARE=Siril+seestar_varstar_siril.py v0.1.5
#FILTER=CV
#DELIM=,
#DATE=JD
#OBSTYPE=CCD
NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
ES UMa,2460325.438000,-0.3450,0.0230,CV,NO,DIFF,ENSEMBLE,na,na,na,na,1,na,seestar_s30pro
"""


class TestAavsoSpecChecker(unittest.TestCase):
    """The checker itself: accepts the spec's examples, rejects the old export."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_spec_examples_pass(self):
        for i, text in enumerate(_SPEC_EXAMPLES, start=1):
            path = self.tmp / f"example{i}.txt"
            path.write_text(text, encoding="utf-8")
            self.assertEqual(aavso_violations(path), [], f"spec example {i}")

    def test_old_export_fails(self):
        path = self.tmp / "aavso.csv"
        path.write_text(_OLD_EXPORT, encoding="utf-8")
        problems = "\n".join(aavso_violations(path))
        self.assertIn("MAG 'MAG' is not a magnitude", problems)  # uncommented header
        self.assertIn("MTYPE 'DIFF'", problems)


class TestExportAavso(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _read_csv_rows(self, path: Path):
        """Data rows keyed by the '#NAME,DATE,...' column-name line."""
        lines = path.read_text(encoding="utf-8").splitlines()
        names = next(l for l in lines if l.startswith("#NAME,"))[1:].split(",")
        return list(csv.DictReader((l for l in lines if not l.startswith("#")),
                                   fieldnames=names))

    def _read_headers(self, path: Path):
        return [l.strip() for l in path.read_text().splitlines()
                if l.startswith("#")]

    def test_valid_rows_written(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.0230",
            "2460325.439000 12.3670 0.0310",
        ])
        out = self.tmp / "out.csv"
        rows = _export(dat, out)
        self.assertEqual(len(rows), 2)
        self.assertTrue(out.exists())

    def test_file_meets_aavso_spec(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 -0.3450 0.0230",
            "2460325.439000 -0.3670 0.0310",
        ])
        out = self.tmp / "aavso.csv"
        _export(dat, out, "V* ES UMa", filt_note="LP_filter_Seestar_S30Pro",
                airmass_map=[(2460325.438, 1.2345), (2460325.439, 1.2351)])
        self.assertEqual(aavso_violations(out), [])

    def test_column_names_line_is_commented(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 -0.3450 0.0230"])
        out = self.tmp / "out.csv"
        _export(dat, out)
        lines = out.read_text().splitlines()
        self.assertIn("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,"
                      "KNAME,KMAG,AMASS,GROUP,CHART,NOTES", lines)
        self.assertFalse(any(l.startswith("NAME,") for l in lines))

    def test_magnitude_is_apparent_ensemble_std(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 -0.3450 0.0230"])
        out = self.tmp / "out.csv"
        _export(dat, out, ensemble_vmag=12.0)
        row = self._read_csv_rows(out)[0]
        self.assertEqual(row["MAG"], "11.655")
        self.assertEqual(row["MERR"], "0.023")
        self.assertEqual((row["MTYPE"], row["TRANS"]), ("STD", "NO"))
        self.assertEqual((row["CNAME"], row["CMAG"]), ("ENSEMBLE", "na"))
        self.assertEqual((row["GROUP"], row["CHART"]), ("na", "APASS DR9"))

    def test_nan_magnitude_excluded(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.0230",
            "2460325.439000 nan 0.0",
        ])
        out = self.tmp / "out.csv"
        rows = _export(dat, out)
        self.assertEqual(len(rows), 1)

    def test_high_merr_excluded(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.0230",
            "2460325.439000 12.3670 0.6000",  # MERR > 0.5
        ])
        out = self.tmp / "out.csv"
        rows = _export(dat, out)
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(rows[0][2], 0.0230, places=4)

    def test_merr_at_threshold_excluded(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.5001",
        ])
        rows = _export(dat, self.tmp / "out.csv")
        self.assertEqual(len(rows), 0)

    def test_merr_at_threshold_included(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.4999",
        ])
        rows = _export(dat, self.tmp / "out.csv")
        self.assertEqual(len(rows), 1)

    def test_empty_dat_returns_empty(self):
        dat = _write_dat(self.tmp / "lc.dat", [])
        rows = _export(dat, self.tmp / "out.csv")
        self.assertEqual(rows, [])
        self.assertFalse((self.tmp / "out.csv").exists())

    def test_missing_obscode_writes_nothing(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 -0.3450 0.0230"])
        out = self.tmp / "out.csv"
        with self.assertRaisesRegex(ValueError, "observer code"):
            _export(dat, out, obscode="  ")
        self.assertFalse(out.exists())

    def test_differential_only_writes_nothing(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 -0.3450 0.0230"])
        out = self.tmp / "out.csv"
        with self.assertRaisesRegex(ValueError, "differential"):
            _export(dat, out, ensemble_vmag=None)
        self.assertFalse(out.exists())

    def test_unknown_filter_writes_nothing(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 -0.3450 0.0230"])
        out = self.tmp / "out.csv"
        with self.assertRaisesRegex(ValueError, "filter"):
            _export(dat, out, filt_code="SRJ")
        self.assertFalse(out.exists())

    def test_csv_headers_present(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        _export(dat, out, obscode="tst01")
        headers = self._read_headers(out)
        self.assertIn("#TYPE=EXTENDED", headers)
        self.assertIn("#OBSCODE=TST01", headers)
        self.assertIn("#DELIM=,", headers)
        self.assertIn("#DATE=JD", headers)
        self.assertIn("#OBSTYPE=CCD", headers)
        self.assertTrue(any("SOFTWARE" in h for h in headers))
        self.assertTrue(any(VERSION in h for h in headers))

    def test_filter_code_in_csv(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        _export(dat, out, filt_code="V")
        self.assertEqual(self._read_csv_rows(out)[0]["FILT"], "V")

    def test_star_name_in_data_rows(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        _export(dat, out)
        rows = self._read_csv_rows(out)
        self.assertEqual(rows[0]["NAME"], "ES UMa")

    def test_vsx_prefix_stripped_from_name(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        _export(dat, out, "V* ES UMa")
        self.assertEqual(self._read_csv_rows(out)[0]["NAME"], "ES UMa")

    def test_filt_note_in_notes_column(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        _export(dat, out, filt_note="LP_filter")
        rows = self._read_csv_rows(out)
        self.assertIn("LP_filter", rows[0]["NOTES"])
        self.assertIn("5 APASS DR9 comparison stars", rows[0]["NOTES"])

    def test_notes_never_contain_the_delimiter(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        _export(dat, out, filt_note="LP, 2 inch")
        self.assertEqual(aavso_violations(out), [])

    def test_unix_line_endings(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        _export(dat, out)
        self.assertNotIn(b"\r", out.read_bytes())

    def test_airmass_written(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 12.3450 0.0230"])
        out = self.tmp / "out.csv"
        _export(dat, out, airmass_map=[(2460325.438, 1.2345)])
        self.assertEqual(self._read_csv_rows(out)[0]["AMASS"], "1.234")

    def test_custom_merr_max(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 12.3450 0.8000",
        ])
        rows = _export(dat, self.tmp / "out.csv", merr_max=1.0)
        self.assertEqual(len(rows), 1)

    def test_every_filter_option_is_an_aavso_code(self):
        for label, (code, _) in FILTER_OPTIONS.items():
            self.assertIn(code, FILTERS, label)


class TestEnsembleZeroPoint(unittest.TestCase):

    def _siril_v_minus_c(self, target_inst: float, comps_inst: list[float]) -> float:
        """V-C as Siril's light_curve computes it (photometry.c, new_light_curve)."""
        flux = sum(10 ** (-0.4 * m) for m in comps_inst) / len(comps_inst)
        return target_inst - (-2.5 * math.log10(flux))

    def test_recovers_the_target_magnitude(self):
        comps_v = [10.2, 11.0, 11.9, 12.8, 13.1]
        target_v, zp = 11.73, -21.4  # instrumental = catalogue + zero point
        vc = self._siril_v_minus_c(target_v + zp, [v + zp for v in comps_v])
        self.assertAlmostEqual(vc + ensemble_zero_point(comps_v), target_v, places=9)
        # The median the pipeline used before is off by ~0.6 mag here.
        self.assertGreater(abs(vc + statistics.median(comps_v) - target_v), 0.5)

    def test_equal_magnitudes(self):
        self.assertAlmostEqual(ensemble_zero_point([12.0, 12.0, 12.0]), 12.0)

    def test_empty(self):
        self.assertIsNone(ensemble_zero_point([]))


class TestSirilLog(unittest.TestCase):
    COMPS = [("Comp1", "1", 11.0), ("Comp1", "11", 12.0), ("Comp1", "2", 13.0)]

    def test_english(self):
        out = [
            "log: Target star identified: ES UMa",
            "log: Star 1 could not be used because it's on the borders or outside",
            "log: star 11 [Comp1] added as a reference star",
            "log: star 2 [Comp1] added as a reference star",
            "log: Using 2 stars to calibrate the light curve",
        ]
        self.assertEqual(siril_loaded_comps(out, self.COMPS), self.COMPS[1:])
        self.assertEqual(siril_refs_used(out), 2)

    def test_french(self):
        out = [
            "log: étoile 1 [Comp1] ajoutée comme étoile de référence",
            "log: Utilisation de 1 étoiles pour étalonner la courbe de lumière",
        ]
        self.assertEqual(siril_loaded_comps(out, self.COMPS), self.COMPS[:1])
        self.assertEqual(siril_refs_used(out), 1)

    def test_single_reference(self):
        out = ["log: Only one reference star was validated, this will not result "
               "in an accurate light curve."]
        self.assertEqual(siril_refs_used(out), 1)

    def test_not_found(self):
        self.assertIsNone(siril_refs_used(["log: something else"]))
        self.assertEqual(siril_loaded_comps(["log: something else"], self.COMPS), [])


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

    def test_parse_comp_csv(self):
        f = _write_comp_csv(self.tmp / "comp.csv", n_comp=2)
        self.assertEqual(_parse_comp_csv(f),
                         [("Comp1", "star1", 12.1), ("Comp1", "star2", 12.2)])


if __name__ == "__main__":
    unittest.main()
