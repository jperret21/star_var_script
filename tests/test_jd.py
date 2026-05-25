"""Tests for JD conversion and light_curve.dat normalization."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import shutil
import tempfile
import unittest
from pathlib import Path

from pipeline import dateobs_to_jd, inject_jd_into_dat
from tests.helpers import make_fits, make_seq, make_dat


class TestDateObsToJD(unittest.TestCase):

    def test_j2000_epoch(self):
        # J2000.0 reference point must map to exactly JD 2451545.0
        jd = dateobs_to_jd("2000-01-01T12:00:00")
        self.assertAlmostEqual(jd, 2451545.0, places=6)

    def test_with_microseconds(self):
        jd = dateobs_to_jd("2024-01-15T22:30:45.500")
        self.assertIsNotNone(jd)
        self.assertGreater(jd, 2451545.0)

    def test_without_microseconds(self):
        jd = dateobs_to_jd("2024-01-15T22:30:45")
        self.assertIsNotNone(jd)

    def test_without_seconds(self):
        jd = dateobs_to_jd("2024-01-15T22:30")
        self.assertIsNotNone(jd)

    def test_one_minute_apart(self):
        jd1 = dateobs_to_jd("2024-01-15T22:30:00")
        jd2 = dateobs_to_jd("2024-01-15T22:31:00")
        self.assertAlmostEqual(jd2 - jd1, 60 / 86400, places=9)

    def test_invalid_returns_none(self):
        self.assertIsNone(dateobs_to_jd("not-a-date"))
        self.assertIsNone(dateobs_to_jd(""))
        self.assertIsNone(dateobs_to_jd("2024-01-15"))  # date only, no T


class TestInjectJD(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    # ── Format A: #JD_UT (+ 0) — frame indices → inject from FITS ────────────

    def test_format_a_injects_real_jd(self):
        stem = "r_light"
        make_fits(self.tmp / f"{stem}_00001.fit",
                  {"DATE-OBS": "2000-01-01T12:00:00", "EXPTIME": "20.0"})
        make_fits(self.tmp / f"{stem}_00002.fit",
                  {"DATE-OBS": "2000-01-01T12:10:00", "EXPTIME": "20.0"})
        seq = make_seq(self.tmp / f"{stem}_.seq", selected=[1, 2], all_count=2)
        dat = make_dat(self.tmp / "light_curve.dat",
                       "#JD_UT (+ 0)",
                       ["1 12.345 0.023", "2 12.367 0.031"])

        result = inject_jd_into_dat(dat, seq, stem, fixlen=5,
                                    proc=self.tmp, exptime_s=20.0)
        self.assertTrue(result)

        lines = [l for l in dat.read_text().splitlines()
                 if l and not l.startswith("#")]
        jd1 = float(lines[0].split()[0])
        jd2 = float(lines[1].split()[0])
        # Both should be valid JD values
        self.assertGreater(jd1, 2_400_000)
        self.assertGreater(jd2, jd1)
        # 10 minutes apart = 10/1440 days
        self.assertAlmostEqual(jd2 - jd1, 10 / 1440, places=4)

    def test_format_a_missing_fits_returns_false(self):
        seq = make_seq(self.tmp / "r_light_.seq", selected=[1], all_count=1)
        dat = make_dat(self.tmp / "light_curve.dat",
                       "#JD_UT (+ 0)", ["1 12.345 0.023"])
        # No FITS files written → jd_map empty → False
        result = inject_jd_into_dat(dat, seq, "r_light", fixlen=5,
                                    proc=self.tmp, exptime_s=20.0)
        self.assertFalse(result)

    # ── Format B: #JD_UT (+ N) — fractional offsets → add N ─────────────────

    def test_format_b_adds_julian0(self):
        julian0 = 2460325
        dat = make_dat(self.tmp / "light_curve.dat",
                       f"#JD_UT (+ {julian0})",
                       ["0.4375 12.345 0.023", "0.4395 12.367 0.031"])
        seq = self.tmp / "r_light_.seq"  # not read for format B
        seq.write_text("")

        result = inject_jd_into_dat(dat, seq, "r_light", fixlen=5,
                                    proc=self.tmp, exptime_s=20.0)
        self.assertTrue(result)

        lines = [l for l in dat.read_text().splitlines()
                 if l and not l.startswith("#")]
        jd1 = float(lines[0].split()[0])
        jd2 = float(lines[1].split()[0])
        self.assertAlmostEqual(jd1, julian0 + 0.4375, places=6)
        self.assertAlmostEqual(jd2, julian0 + 0.4395, places=6)
        self.assertGreater(jd2, jd1)

    # ── Format C: #JD_UT — already absolute JD → no-op ───────────────────────

    def test_format_c_already_absolute_is_noop(self):
        original = "#JD_UT\n2460325.438000 12.345 0.023\n"
        dat = self.tmp / "light_curve.dat"
        dat.write_text(original, encoding="utf-8")
        seq = self.tmp / "r_light_.seq"
        seq.write_text("")

        result = inject_jd_into_dat(dat, seq, "r_light", fixlen=5,
                                    proc=self.tmp, exptime_s=20.0)
        self.assertFalse(result)
        self.assertEqual(dat.read_text(), original)

    # ── Half-exposure offset ──────────────────────────────────────────────────

    def test_format_a_applies_half_exptime(self):
        stem = "r_light"
        make_fits(self.tmp / f"{stem}_00001.fit",
                  {"DATE-OBS": "2000-01-01T12:00:00", "EXPTIME": "20.0"})
        seq = make_seq(self.tmp / f"{stem}_.seq", selected=[1], all_count=1)
        dat = make_dat(self.tmp / "light_curve.dat",
                       "#JD_UT (+ 0)", ["1 12.345 0.023"])

        inject_jd_into_dat(dat, seq, stem, fixlen=5,
                           proc=self.tmp, exptime_s=20.0)
        line = [l for l in dat.read_text().splitlines()
                if l and not l.startswith("#")][0]
        jd = float(line.split()[0])
        jd_start = dateobs_to_jd("2000-01-01T12:00:00")
        # JD should be start + 10s (half of 20s exptime).
        # places=5 (≈ 0.9s precision) accounts for the 6-decimal rounding in the .dat file.
        self.assertAlmostEqual(jd - jd_start, 10 / 86400, places=5)


if __name__ == "__main__":
    unittest.main()
