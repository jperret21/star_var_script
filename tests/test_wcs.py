"""Tests for sky_to_pixel, stars_in_safe_circle, stars_in_frame and plate_scale_arcsec."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import math
import unittest

from pipeline import plate_scale_arcsec, sky_to_pixel, stars_in_frame, stars_in_safe_circle


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

    SEESTAR_CDELT = 2.9e-3 / 160 * (180 / math.pi)  # 1.0385e-3 deg/px = 3.74 arcsec/px

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


class TestSkyToPixelRotated(unittest.TestCase):
    """sky_to_pixel with a real rotated WCS (PC matrix ≠ identity).

    Values taken from an actual Seestar plate-solved r_light frame (~42° rotation).
    This exercises the rotated-field code path where CDELT × PC matrix ≠ simple
    diagonal CD matrix.
    """

    HDR = {
        "NAXIS1": 2160.0, "NAXIS2": 3840.0,
        "CRVAL1": 149.246533, "CRVAL2": 69.072210,
        "CRPIX1": 1080.5,    "CRPIX2": 1920.5,
        "CDELT1": -0.00102032, "CDELT2": 0.00102019,
        "PC1_1": -0.736069,  "PC1_2":  0.676906,
        "PC2_1":  0.676688,  "PC2_2":  0.736270,
    }

    def _in_frame(self, ra, dec, margin=35):
        naxis1 = int(self.HDR["NAXIS1"]); naxis2 = int(self.HDR["NAXIS2"])
        px, py = sky_to_pixel(ra, dec, self.HDR)
        if px is None:
            return False
        tdx = round(px - 0.5)
        tdy = round(naxis2 - py + 0.5)
        return (margin < tdx < naxis1 - margin and
                margin < tdy < naxis2 - margin)

    def test_crval_maps_to_crpix(self):
        px, py = sky_to_pixel(self.HDR["CRVAL1"], self.HDR["CRVAL2"], self.HDR)
        self.assertAlmostEqual(px, self.HDR["CRPIX1"], places=3)
        self.assertAlmostEqual(py, self.HDR["CRPIX2"], places=3)

    def test_edge_star_out_of_frame(self):
        # NSVS J0941171+685517: 3.9° from field centre in RA. Under the correct
        # WCS it projects to disp-x ≈ -8 (off the left edge), so it is NOT in the
        # frame. (Cross-check: astropy wcs_world2pix agrees to <0.01 px, and
        # test_nsvs_edge_star_excluded excludes the same star from the safe circle.)
        self.assertFalse(self._in_frame(145.3213, 68.9214))

    def test_apass_comps_near_edge_star_excluded(self):
        # APASS stars found within 1.5° of the edge star project BELOW the frame
        # (disp-y > NAXIS2) because the field rotation sends them out of bounds.
        out_of_frame_samples = [
            (143.7988, 67.6621),
            (143.0804, 68.0335),
            (143.1283, 67.7389),
        ]
        for ra, dec in out_of_frame_samples:
            self.assertFalse(self._in_frame(ra, dec),
                             msg=f"Expected ({ra},{dec}) to be outside frame")

    def test_field_centre_comps_in_frame(self):
        # APASS stars around CRVAL are guaranteed to be in frame.
        centre_samples = [
            (149.7710, 68.7973),
            (148.9206, 69.0633),
            (149.2537, 68.9019),
        ]
        for ra, dec in centre_samples:
            self.assertTrue(self._in_frame(ra, dec),
                            msg=f"Expected ({ra},{dec}) to be inside frame")


class TestStarsInSafeCircle(unittest.TestCase):
    """stars_in_safe_circle — alt-az inscribed circle filter."""

    # Real Seestar rotated WCS (same header as TestSkyToPixelRotated)
    HDR    = TestSkyToPixelRotated.HDR
    W, H   = int(HDR["NAXIS1"]), int(HDR["NAXIS2"])   # 2160 × 3840
    MARGIN = 50

    def _star(self, ra, dec, name="S"):
        return {"name": name, "ra": ra, "dec": dec, "mag": 12.0}

    # ── Basic invariants ──────────────────────────────────────────────────────

    def test_empty_input_returns_empty(self):
        self.assertEqual(stars_in_safe_circle([], self.HDR, self.W, self.H), [])

    def test_no_wcs_returns_all(self):
        # Without WCS we cannot project → return stars unchanged (conservative)
        s = self._star(self.HDR["CRVAL1"], self.HDR["CRVAL2"])
        result = stars_in_safe_circle([s], {}, self.W, self.H)
        self.assertEqual(len(result), 1)

    def test_crval_star_included(self):
        # Star at field centre is always inside the safe circle
        s = self._star(self.HDR["CRVAL1"], self.HDR["CRVAL2"])
        result = stars_in_safe_circle([s], self.HDR, self.W, self.H, self.MARGIN)
        self.assertEqual(len(result), 1)

    # ── Real-world cases ──────────────────────────────────────────────────────

    def test_nsvs_edge_star_excluded(self):
        # NSVS J0941171+685517: ~1387 px from centre > safe_r=1030 px → excluded
        s = self._star(145.3213, 68.9214, "NSVS")
        result = stars_in_safe_circle([s], self.HDR, self.W, self.H, self.MARGIN)
        self.assertEqual(len(result), 0)

    def test_es_uma_included(self):
        # ES UMa: ~264 px from centre << safe_r=1030 px → included
        s = self._star(148.619236, 69.222850, "ES UMa")
        result = stars_in_safe_circle([s], self.HDR, self.W, self.H, self.MARGIN)
        self.assertEqual(len(result), 1)

    def test_apass_field_centre_comps_included(self):
        # APASS stars near CRVAL (confirmed in-frame by earlier tests)
        centre_samples = [
            (149.7710, 68.7973),
            (148.9206, 69.0633),
            (149.2537, 68.9019),
        ]
        stars = [self._star(ra, dec, f"S{i}") for i, (ra, dec) in enumerate(centre_samples)]
        result = stars_in_safe_circle(stars, self.HDR, self.W, self.H, self.MARGIN)
        self.assertEqual(len(result), 3)

    def test_mixed_in_and_out(self):
        centre = self._star(self.HDR["CRVAL1"], self.HDR["CRVAL2"], "centre")
        edge   = self._star(145.3213, 68.9214, "NSVS")
        result = stars_in_safe_circle([centre, edge], self.HDR, self.W, self.H, self.MARGIN)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["name"], "centre")

    # ── Radius boundary ───────────────────────────────────────────────────────

    def test_safe_radius_is_min_half_minus_margin(self):
        # Verify computed safe_r = min(W,H)/2 - margin = 1080 - 50 = 1030 px
        safe_r = min(self.W, self.H) / 2.0 - self.MARGIN
        self.assertAlmostEqual(safe_r, 1030.0, places=1)

    def test_zero_margin_gives_full_inscribed_circle(self):
        # With margin=0, NSVS at 1387 px should still be outside (>1080)
        s = self._star(145.3213, 68.9214, "NSVS")
        result = stars_in_safe_circle([s], self.HDR, self.W, self.H, margin=0)
        self.assertEqual(len(result), 0)



class TestPlateScale(unittest.TestCase):

    def test_cdelt_header(self):
        hdr = _cdelt_hdr(291.4, 41.9, 917.5, 1126.0, 1.0204e-3)
        self.assertAlmostEqual(plate_scale_arcsec(hdr), 3.673, places=3)

    def test_cd_header(self):
        hdr = _cd_hdr(291.4, 41.9, 461.5, 562.0, 2.0413e-3)
        self.assertAlmostEqual(plate_scale_arcsec(hdr), 7.349, places=3)

    def test_rotation_does_not_change_scale(self):
        c, s = math.cos(math.radians(3.3)), math.sin(math.radians(3.3))
        hdr = _cdelt_hdr(291.4, 41.9, 917.5, 1126.0, 1.0204e-3)
        hdr.update({"PC1_1": c, "PC1_2": -s, "PC2_1": -s, "PC2_2": -c})
        self.assertAlmostEqual(plate_scale_arcsec(hdr), 3.673, places=3)

    def test_no_wcs(self):
        self.assertIsNone(plate_scale_arcsec({}))



class TestStarsInFrame(unittest.TestCase):
    """stars_in_frame — equatorial mode, whole frame minus an edge margin."""

    # Real Seestar S30 Pro registered frame (-framing=min), RR Lyr session shot
    # in equatorial mode, 2026-09-25.
    HDR = {
        "NAXIS1": 1834.0, "NAXIS2": 2251.0,
        "CRVAL1": 291.444118262612, "CRVAL2": 41.9003498874388,
        "CRPIX1": 917.5, "CRPIX2": 1126.0,
        "CDELT1": -0.00102037074573893, "CDELT2": 0.00102053780890616,
        "PC1_1": 0.998319905047237, "PC1_2": -0.0579427923600136,
        "PC2_1": -0.0580727356630557, "PC2_2": -0.998312354612828,
    }
    W, H = int(HDR["NAXIS1"]), int(HDR["NAXIS2"])
    RR_LYR = {"name": "RR Lyr", "ra": 291.36630, "dec": 42.78436, "mag": 7.2}

    def test_no_wcs_returns_all(self):
        self.assertEqual(len(stars_in_frame([self.RR_LYR], {}, self.W, self.H)), 1)

    def test_rr_lyr_outside_safe_circle_but_inside_frame(self):
        # ~869 px from the centre, 257 px from the frame edge: rejected by the
        # alt-az safe circle (r = 867 px), kept in equatorial mode.
        self.assertEqual(stars_in_safe_circle([self.RR_LYR], self.HDR, self.W, self.H), [])
        self.assertEqual(stars_in_frame([self.RR_LYR], self.HDR, self.W, self.H), [self.RR_LYR])

    def test_star_off_the_frame_excluded(self):
        # 000-BCG-953 projects to display x ≈ -26, left of the frame.
        s = {"name": "000-BCG-953", "ra": 292.67158333, "dec": 43.03772222}
        self.assertEqual(stars_in_frame([s], self.HDR, self.W, self.H), [])

    def test_margin_excludes_star_near_edge(self):
        # RR Lyr is 257 px from the nearest edge.
        self.assertEqual(stars_in_frame([self.RR_LYR], self.HDR, self.W, self.H, margin=250),
                         [self.RR_LYR])
        self.assertEqual(stars_in_frame([self.RR_LYR], self.HDR, self.W, self.H, margin=260), [])


if __name__ == "__main__":
    unittest.main()
