"""Physical constraints of the public canopy shortwave closure (no reference retuning)."""
import math
import unittest

import numpy as np

from base import canopy_closure as cc
from base import node_radiosity as nr


CASES = (("hom", cc.HOMOGENEOUS), ("row", (0.5, 0.5, 1.0)),
         ("crown", (0.02, 3.0, 1.5)))


def hemisphere(kind, shape, optics, *, diffuse=0.0, sunlit=True, order=48, azimuths=96):
    mu, weights = np.polynomial.legendre.leggauss(order)
    integral = 0.0
    for m, w in zip((mu + 1.0) / 2.0, weights / 2.0):
        z = math.degrees(math.acos(m))
        values = [cc.reflectance(kind, 2.0, 0.1, shape, *optics, z, a, 40.0, 130.0,
                                 diffuse=diffuse, sunlit=sunlit)
                  for a in (np.arange(azimuths) + 0.5) * 360.0 / azimuths]
        integral += 2.0 * w * m * np.mean(values)
    return float(integral)


class OpticalConservationTests(unittest.TestCase):
    def test_leaf_transmission_cannot_backscatter_at_first_order(self):
        for kind, shape in CASES:
            total = cc.reflectance(kind, 2., .1, shape, 0., 1., 0., 40., 130., 40., 130.)
            multiple = nr.shortwave_multiple(cc.NODES[kind], 2., cc._shape(kind, shape), 0., 1., 0.,
                                             40., 130., 40., 130., 0., 1., 0.)
            self.assertAlmostEqual(total - multiple, 0., delta=1.e-12)

    def test_row_grazing_view_has_finite_nonnegative_reflectance(self):
        value = cc.reflectance("row", 2., .1, (.5, .5, 1.), .45, .4, .25,
                               89.999, 45., 40., 130.)
        self.assertTrue(math.isfinite(value))
        self.assertGreaterEqual(value, 0.)

    def test_zero_terrain_density_is_the_plane_canopy(self):
        from base.optical_terrain_standard import reflectance
        for kind, shape in CASES:
            settings = dict(vegetation_model=kind, lai=2., hspot=.1, terrain_density=0.,
                            stand_density=.02, height=3., width=1.5)
            expected = cc.reflectance(kind, 2., .1, shape, .45, .4, .25, 70., 280., 40., 130.)
            for density, tolerance in ((1.e-10, 1.e-4), (0., 1.e-12)):
                actual = reflectance(dict(settings, terrain_density=density),
                                     40., 130., 70., 280., .45, .4, .25)
                self.assertAlmostEqual(actual, expected, delta=tolerance)

    def test_hotspot_width_changes_directional_reflectance(self):
        # At the exact hotspot the joint gap is P(s) for any hspot; the width shows next to it.
        for kind, shape in CASES:
            narrow = cc.reflectance(kind, 2., .05, shape, .45, .40, .25, 45., 130., 40., 130.)
            broad = cc.reflectance(kind, 2., .5, shape, .45, .40, .25, 45., 130., 40., 130.)
            self.assertGreater(abs(narrow - broad), 1.e-5)

    # Known limitations of the shared four-component optical first order (TiRT theory doc 4.7):
    # the standard sunlit-leaf two-layer term and the crown leaf term are not sun/view symmetric,
    # and the directional first order escapes more than the node budget assumes.
    def _reciprocity(self, cases):
        for kind, shape, level in cases:
            with self.subTest(kind=kind, level=level):
                a = cc.reflectance(kind, 2.0, 0.1, shape, .45, .40, .25,
                                   70., 280., 40., 130., row_azimuth=17., level=level)
                b = cc.reflectance(kind, 2.0, 0.1, shape, .45, .40, .25,
                                   40., 130., 70., 280., row_azimuth=17., level=level)
                self.assertAlmostEqual(a, b, delta=1.e-12)

    def test_sun_view_reciprocity(self):
        self._reciprocity([(kind, shape, "upgrade") for kind, shape in CASES if kind != "crown"])

    @unittest.expectedFailure
    def test_sun_view_reciprocity_known_failures(self):
        self._reciprocity([(kind, shape, "standard") for kind, shape in CASES] + [("crown", CASES[2][1], "upgrade")])

    def test_nonabsorbing_sky_reflects_unit_incident_flux(self):
        for kind, shape in CASES:
            with self.subTest(kind=kind):
                self.assertAlmostEqual(hemisphere(kind, shape, (.5, .5, 1.), diffuse=1.0), 1.0, delta=5.e-4)

    @unittest.expectedFailure
    def test_nonabsorbing_beam_reflects_unit_incident_flux(self):
        for kind, shape in CASES:
            with self.subTest(kind=kind):
                self.assertAlmostEqual(hemisphere(kind, shape, (.5, .5, 1.), diffuse=0.0), 1.0, delta=5.e-4)

    def _budget(self, diffuses):
        for kind, shape in CASES:
            nodes, beam = nr.canopy_nodes(cc.NODES[kind], 2.0, shape, 40., 130.)
            for diffuse in diffuses:
                with self.subTest(kind=kind, diffuse=diffuse):
                    arrivals = ((1. - diffuse) * nr.from_above(nodes, beam)
                                + diffuse * nr.sky_arrivals(nodes))
                    first, multiple = nr.absorbed(nodes, .25, .45, .40, arrivals)
                    # Closed opaque-ground energy equation, not another BRF calculation.
                    expected = 1.0 - float((first + multiple).sum())
                    actual = hemisphere(kind, shape, (.45, .40, .25), diffuse=diffuse)
                    self.assertAlmostEqual(actual, expected, delta=2.e-5)

    def test_sky_hemisphere_matches_independent_node_absorption_budget(self):
        self._budget((1.0,))

    @unittest.expectedFailure
    def test_beam_hemisphere_matches_independent_node_absorption_budget(self):
        self._budget((0.0, .35))

    def test_absorbing_and_shadow_limits(self):
        for kind, shape in CASES:
            for diffuse in (0., .5, 1.):
                with self.subTest(kind=kind, diffuse=diffuse):
                    self.assertEqual(cc.reflectance(kind, 2., .1, shape, 0., 0., 0.,
                                                     25., 90., 40., 130., diffuse=diffuse), 0.)
                    shadow = cc.reflectance(kind, 2., .1, shape, .45, .40, .25,
                                            25., 90., 40., 130., diffuse=diffuse, sunlit=False)
                    sky = cc.reflectance(kind, 2., .1, shape, .45, .40, .25,
                                         25., 90., 40., 130., diffuse=1.)
                    self.assertAlmostEqual(shadow, diffuse * sky, delta=1.e-12)

    def test_zero_canopy_returns_bare_ground(self):
        for kind, shape in CASES:
            self.assertAlmostEqual(cc.reflectance(kind, 0., .1, shape, .45, .4, .25,
                                                  70., 280., 40., 130.), .25)
        self.assertAlmostEqual(cc.reflectance("crown", 2., .1, (0., 3., 1.5), .45, .4, .25,
                                              70., 280., 40., 130.), .25)

    def test_zero_crown_density_has_bare_soil_energy_budget(self):
        # Conservative disconnected canopy states must not enter a singular solve.
        shape = (0., 3., 1.5)
        self.assertEqual(cc.beam_albedo("crown", 2., shape, .5, .5, .25, 40., 130.), .25)
        self.assertEqual(cc.sky_albedo("crown", 2., shape, .5, .5, .25), .25)


if __name__ == "__main__":
    unittest.main()
