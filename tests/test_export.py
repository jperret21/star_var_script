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
    _read_lc_dat,
    argos_selection_stars,
    ensemble_zero_point,
    export_aavso,
    load_argos_selection,
    siril_loaded_comps,
    siril_refs_used,
    truncate_comp_csv,
    write_nina_csv,
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

    def test_check_star_and_vsp_chart(self):
        dat = _write_dat(self.tmp / "lc.dat", [
            "2460325.438000 -0.3450 0.0230",
            "2460325.439000 -0.3670 0.0310",
        ])
        out = self.tmp / "aavso.csv"
        _export(dat, out, comp_source="AAVSO VSP", chart="X42585ESI",
                check=("000-BJV-171", {2460325.438: 12.4567}))
        rows = self._read_csv_rows(out)
        self.assertEqual((rows[0]["KNAME"], rows[0]["KMAG"]), ("000-BJV-171", "12.457"))
        self.assertEqual((rows[1]["KNAME"], rows[1]["KMAG"]), ("na", "na"))  # no check point
        self.assertEqual(rows[0]["CHART"], "X42585ESI")
        self.assertIn("5 AAVSO VSP comparison stars", rows[0]["NOTES"])
        self.assertIn("check star 000-BJV-171", rows[0]["NOTES"])
        self.assertEqual(aavso_violations(out), [])

    def test_unknown_chart_is_na(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 -0.3450 0.0230"])
        out = self.tmp / "aavso.csv"
        _export(dat, out, chart=None)
        self.assertEqual(self._read_csv_rows(out)[0]["CHART"], "na")


def _manifest(comps, checks=(), target=(148.619236, 69.222850)):
    """An Argos photometry_selection.json, as argos TargetSet writes it."""
    def rec(auid, ra, dec, v=None, chart="X42585ESI"):
        return {"auid": auid, "name": auid, "ra_deg_j2000": ra, "dec_deg_j2000": dec,
                "catalogue_source": "vsp_auto",
                "catalogue_magnitudes": {} if v is None else {"V": v, "B": v + 0.6},
                "catalogue_chart_id": chart}
    return {"schema": 1, "object_name": "ES UMa",
            "targets": [{"name": "ES UMa", "auid": "000-BCK-001",
                         "ra_deg_j2000": target[0], "dec_deg_j2000": target[1]}],
            "comparison_stars": [rec(*c) for c in comps],
            "check_stars": [rec(*c) for c in checks]}


class TestArgosSelection(unittest.TestCase):
    STAR = {"name": "V* ES UMa", "ra": 148.619236, "dec": 69.222850, "mag": 12.0}

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_comps_check_and_chart(self):
        sel = _manifest(comps=[("000-AAA-001", 148.7, 69.3, 11.4),
                               ("000-AAA-002", 148.5, 69.1, 12.1),
                               ("000-AAA-003", 148.6, 69.2, None),       # no V
                               ("000-AAA-004", 160.0, 69.2, 12.0)],      # off frame
                        checks=[("000-AAA-005", 148.8, 69.2, 12.6)])
        comps, check, chart = argos_selection_stars(
            sel, self.STAR, in_frame=lambda ra, dec: ra < 150)
        self.assertEqual([c[0] for c in comps], ["000-AAA-001", "000-AAA-002"])
        self.assertEqual(comps[0][3], 11.4)
        self.assertEqual(check, ("000-AAA-005", 148.8, 69.2))
        self.assertEqual(chart, "X42585ESI")

    def test_mixed_charts_leave_no_chart(self):
        sel = _manifest(comps=[("000-AAA-001", 148.7, 69.3, 11.4, "X1"),
                               ("000-AAA-002", 148.5, 69.1, 12.1, "X2")])
        _, check, chart = argos_selection_stars(sel, self.STAR, in_frame=lambda *_: True)
        self.assertIsNone(chart)
        self.assertIsNone(check)

    def test_other_target_is_refused(self):
        sel = _manifest(comps=[("000-AAA-001", 148.7, 69.3, 11.4),
                               ("000-AAA-002", 148.5, 69.1, 12.1)],
                        target=(150.0, 69.0))
        with self.assertRaisesRegex(ValueError, "is for ES UMa"):
            argos_selection_stars(sel, self.STAR, in_frame=lambda *_: True)

    def test_too_few_comps_is_refused(self):
        sel = _manifest(comps=[("000-AAA-001", 148.7, 69.3, 11.4)])
        with self.assertRaisesRegex(ValueError, "only 1"):
            argos_selection_stars(sel, self.STAR, in_frame=lambda *_: True)

    def test_load_and_nina_round_trip(self):
        path = self.tmp / "photometry_selection.json"
        path.write_text(__import__("json").dumps(_manifest(
            comps=[("000-AAA-001", 148.7, 69.3, 11.4), ("000-AAA-002", 148.5, 69.1, 12.1)])))
        sel = load_argos_selection(path)
        comps, _, _ = argos_selection_stars(sel, self.STAR, in_frame=lambda *_: True)
        nina = self.tmp / "comp_stars.csv"
        write_nina_csv(nina, ("V* ES UMa", 148.619236, 69.22285, 12.0), comps)
        self.assertTrue(nina.read_text().splitlines()[1].startswith("Target,V* ES UMa,"))
        self.assertEqual(_parse_comp_csv(nina),
                         [("Comp2", "000-AAA-001", 11.4), ("Comp2", "000-AAA-002", 12.1)])
        self.assertIsNone(load_argos_selection(self.tmp / "missing.json"))

    def test_siril_log_names_comp2_stars(self):
        comps = [("Comp2", "000-AAA-001", 11.4), ("Comp2", "000-AAA-002", 12.1)]
        out = ["log: star 000-AAA-002 [Comp2] added as a reference star"]
        self.assertEqual(siril_loaded_comps(out, comps), comps[1:])

    def test_read_lc_dat(self):
        dat = _write_dat(self.tmp / "lc.dat", ["2460325.438000 -0.3450 0.0230",
                                               "2460325.439000 nan 9.999"])
        rows = _read_lc_dat(dat)
        self.assertEqual(rows[0], (2460325.438, -0.345, 0.023))
        self.assertTrue(math.isnan(rows[1][1]))


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
