"""Tests for sky_to_pixel and stars_in_frame."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import math
import unittest

from pipeline import sky_to_pixel, stars_in_frame


def _cdelt_hdr(crval1, crval2, crpix1, crpix2, cdelt_deg):
    """Minimal WCS header using CDELT + PC matrix (identity rotation)."""
    return {
        "CRVAL1": crval1, "CRVAL2": crval2,
        "CRPIX1": crpix1, "CRPIX2": crpix2,
        "CDELT1": -cdelt_deg, "CDELT2": cdelt_deg,
        "PC1_1": 1.0, "PC1_2": 0.0,
        "PC2_1": 0.0, "PC2_2": 1.0,
    }


def _cd_hdr(crval1, crval2, crpix1, crpix2, cdelt_deg):
    """Minimal WCS header using CD matrix (identity rotation, no skew)."""
    return {
        "CRVAL1": crval1, "CRVAL2": crval2,
        "CRPIX1": crpix1, "CRPIX2": crpix2,
        "CD1_1": -cdelt_deg, "CD1_2": 0.0,
        "CD2_1": 0.0,        "CD2_2": cdelt_deg,
    }


class TestSkyToPixel(unittest.TestCase):

    SEESTAR_CDELT = 2.9e-3 / 160 * (180 / math.pi)  # ~1.035 arcsec/pix in degrees

    # ── Invariant: projecting CRVAL gives back CRPIX ─────────────────────────

    def test_center_point_cdelt(self):
        hdr = _cdelt_hdr(148.888, 69.065, 1080.0, 1920.0, self.SEESTAR_CDELT)
        px, py = sky_to_pixel(148.888, 69.065, hdr)
        self.assertAlmostEqual(px, 1080.0, places=6)
        self.assertAlmostEqual(py, 1920.0, places=6)

    def test_center_point_cd(self):
        hdr = _cd_hdr(148.888, 69.065, 1080.0, 1920.0, self.SEESTAR_CDELT)
        px, py = sky_to_pixel(148.888, 69.065, hdr)
        self.assertAlmostEqual(px, 1080.0, places=6)
        self.assertAlmostEqual(py, 1920.0, places=6)

    def test_center_equatorial(self):
        # Field centered on the equator
        hdr = _cdelt_hdr(180.0, 0.0, 1080.0, 1920.0, self.SEESTAR_CDELT)
        px, py = sky_to_pixel(180.0, 0.0, hdr)
        self.assertAlmostEqual(px, 1080.0, places=6)
        self.assertAlmostEqual(py, 1920.0, places=6)

    # ── Offset: 1 pixel in RA direction ──────────────────────────────────────

    def test_one_pixel_offset(self):
        cdelt = self.SEESTAR_CDELT
        hdr = _cdelt_hdr(148.888, 69.065, 1080.0, 1920.0, cdelt)
        # delta_ra chosen so that sky_to_pixel returns exactly 1 pixel offset in x
        delta_ra = cdelt / math.cos(math.radians(69.065))
        px, py = sky_to_pixel(148.888 + delta_ra, 69.065, hdr)
        # Verify the offset magnitude is 1 pixel and Dec axis is unchanged
        self.assertAlmostEqual(abs(px - 1080.0), 1.0, places=1)
        self.assertAlmostEqual(py, 1920.0, places=1)

    # ── Missing keywords → (None, None) ──────────────────────────────────────

    def test_empty_header(self):
        px, py = sky_to_pixel(148.888, 69.065, {})
        self.assertIsNone(px)
        self.assertIsNone(py)

    def test_missing_crval(self):
        hdr = {"CRPIX1": 1080, "CRPIX2": 1920,
               "CDELT1": -self.SEESTAR_CDELT, "CDELT2": self.SEESTAR_CDELT}
        px, py = sky_to_pixel(148.888, 69.065, hdr)
        self.assertIsNone(px)
        self.assertIsNone(py)

    def test_singular_cd_matrix(self):
        hdr = {
            "CRVAL1": 148.888, "CRVAL2": 69.065,
            "CRPIX1": 1080.0,  "CRPIX2": 1920.0,
            "CD1_1": 0.0, "CD1_2": 0.0,
            "CD2_1": 0.0, "CD2_2": 0.0,
        }
        px, py = sky_to_pixel(148.888, 69.065, hdr)
        self.assertIsNone(px)
        self.assertIsNone(py)

    # ── Symmetry: star at same RA but ±1° in Dec ─────────────────────────────

    def test_dec_symmetry(self):
        hdr = _cdelt_hdr(180.0, 0.0, 1080.0, 1920.0, self.SEESTAR_CDELT)
        _, py_plus  = sky_to_pixel(180.0,  1.0, hdr)
        _, py_minus = sky_to_pixel(180.0, -1.0, hdr)
        # Symmetric around CRPIX2
        self.assertAlmostEqual(py_plus + py_minus, 2 * 1920.0, places=3)


class TestStarsInFrame(unittest.TestCase):

    CDELT = TestSkyToPixel.SEESTAR_CDELT  # ~1.035 arcsec/px in degrees
    CRVAL1 = 148.888
    CRVAL2 = 69.065
    CRPIX1 = 1080.0
    CRPIX2 = 1920.0
    W, H = 2160, 3840  # typical Seestar frame

    def _hdr(self):
        return _cdelt_hdr(self.CRVAL1, self.CRVAL2,
                          self.CRPIX1, self.CRPIX2, self.CDELT)

    def _star(self, ra, dec, name="S"):
        return {"name": name, "ra": ra, "dec": dec, "mag": 12.0}

    def test_center_star_included(self):
        s = self._star(self.CRVAL1, self.CRVAL2)
        result = stars_in_frame([s], self._hdr(), self.W, self.H)
        self.assertEqual(len(result), 1)

    def test_far_out_of_frame_excluded(self):
        # 5 degrees away — well outside a 2.2×3.9° field
        s = self._star(self.CRVAL1 + 5.0, self.CRVAL2)
        result = stars_in_frame([s], self._hdr(), self.W, self.H)
        self.assertEqual(len(result), 0)

    def test_margin_excludes_edge_star(self):
        # Place a star exactly on pixel (margin-1, CRPIX2) → should be excluded
        margin = 50
        # One pixel inside the margin from the left edge
        cdelt = self.CDELT
        delta_ra = (margin - 1 - self.CRPIX1) * (-cdelt) / math.cos(math.radians(self.CRVAL2))
        s = self._star(self.CRVAL1 + delta_ra, self.CRVAL2)
        result = stars_in_frame([s], self._hdr(), self.W, self.H, margin=margin)
        self.assertEqual(len(result), 0)

    def test_empty_input_returns_empty(self):
        result = stars_in_frame([], self._hdr(), self.W, self.H)
        self.assertEqual(result, [])

    def test_no_wcs_skips_all(self):
        # Without WCS keywords sky_to_pixel returns (None, None) → all excluded
        s = self._star(self.CRVAL1, self.CRVAL2)
        result = stars_in_frame([s], {}, self.W, self.H)
        self.assertEqual(result, [])

    def test_mixed_in_and_out(self):
        center = self._star(self.CRVAL1, self.CRVAL2, "center")
        far    = self._star(self.CRVAL1 + 5.0, self.CRVAL2, "far")
        result = stars_in_frame([center, far], self._hdr(), self.W, self.H)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["name"], "center")


if __name__ == "__main__":
    unittest.main()
