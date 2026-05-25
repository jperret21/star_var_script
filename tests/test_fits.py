"""Tests for FITS header I/O: read, strip, restore."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import shutil
import tempfile
import unittest
from pathlib import Path

from pipeline import read_fits_header, _fits_strip_keyword, _fits_restore_keyword
from tests.helpers import make_fits


class TestReadFitsHeader(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_reads_string_keyword(self):
        f = make_fits(self.tmp / "test.fit", {"BAYERPAT": "GRBG"})
        hdr = read_fits_header(f)
        self.assertEqual(hdr.get("BAYERPAT"), "GRBG")

    def test_reads_numeric_keyword(self):
        f = make_fits(self.tmp / "test.fit", {"EXPTIME": 20.0})
        hdr = read_fits_header(f)
        self.assertAlmostEqual(float(hdr["EXPTIME"]), 20.0)

    def test_reads_date_obs_as_string(self):
        f = make_fits(self.tmp / "test.fit",
                      {"DATE-OBS": "2024-01-15T22:30:45"})
        hdr = read_fits_header(f)
        self.assertEqual(hdr.get("DATE-OBS"), "2024-01-15T22:30:45")

    def test_reads_multiple_keywords(self):
        f = make_fits(self.tmp / "test.fit",
                      {"BAYERPAT": "GRBG", "EXPTIME": 20.0,
                       "DATE-OBS": "2024-01-15T22:30:00"})
        hdr = read_fits_header(f)
        self.assertEqual(hdr.get("BAYERPAT"), "GRBG")
        self.assertAlmostEqual(float(hdr["EXPTIME"]), 20.0)
        self.assertEqual(hdr.get("DATE-OBS"), "2024-01-15T22:30:00")

    def test_missing_file_returns_empty(self):
        hdr = read_fits_header(self.tmp / "nonexistent.fit")
        self.assertEqual(hdr, {})


class TestFitsStripRestore(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_strip_returns_offset_and_value(self):
        f = make_fits(self.tmp / "test.fit", {"BAYERPAT": "GRBG"})
        result = _fits_strip_keyword(f, "BAYERPAT")
        self.assertIsNotNone(result)
        offset, value = result
        self.assertIsInstance(offset, int)
        self.assertEqual(value, "GRBG")

    def test_strip_removes_keyword(self):
        f = make_fits(self.tmp / "test.fit", {"BAYERPAT": "GRBG"})
        _fits_strip_keyword(f, "BAYERPAT")
        hdr = read_fits_header(f)
        self.assertNotIn("BAYERPAT", hdr)

    def test_strip_nonexistent_returns_none(self):
        f = make_fits(self.tmp / "test.fit", {"EXPTIME": 20.0})
        result = _fits_strip_keyword(f, "BAYERPAT")
        self.assertIsNone(result)

    def test_restore_puts_keyword_back(self):
        f = make_fits(self.tmp / "test.fit", {"BAYERPAT": "GRBG"})
        offset, value = _fits_strip_keyword(f, "BAYERPAT")
        # Keyword is gone
        self.assertNotIn("BAYERPAT", read_fits_header(f))
        # Restore
        ok = _fits_restore_keyword(f, "BAYERPAT", offset, value)
        self.assertTrue(ok)
        hdr = read_fits_header(f)
        self.assertEqual(hdr.get("BAYERPAT"), "GRBG")

    def test_roundtrip_preserves_other_keywords(self):
        f = make_fits(self.tmp / "test.fit",
                      {"BAYERPAT": "GRBG", "EXPTIME": 20.0})
        offset, value = _fits_strip_keyword(f, "BAYERPAT")
        _fits_restore_keyword(f, "BAYERPAT", offset, value)
        hdr = read_fits_header(f)
        self.assertEqual(hdr.get("BAYERPAT"), "GRBG")
        self.assertAlmostEqual(float(hdr["EXPTIME"]), 20.0)

    def test_strip_does_not_affect_other_keywords(self):
        f = make_fits(self.tmp / "test.fit",
                      {"BAYERPAT": "GRBG", "DATE-OBS": "2024-01-15T22:30:00"})
        _fits_strip_keyword(f, "BAYERPAT")
        hdr = read_fits_header(f)
        self.assertNotIn("BAYERPAT", hdr)
        self.assertEqual(hdr.get("DATE-OBS"), "2024-01-15T22:30:00")


if __name__ == "__main__":
    unittest.main()
