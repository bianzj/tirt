"""Urban optical energy, shadows and whole-crown source-profile limits."""
import unittest
from unittest.mock import patch

import numpy as np

import run
from base.optical_crown import crown_reflectance, transport, geometry
from base.optical_crown_profile import source_profile
from base.optical_urban import urban_reflectance, wall_exchange


CITY = dict(urban_building_length=10., urban_building_width=10., urban_building_height=20.,
            urban_density=.0015, lai=0., urban_vegetation_fraction=0.)
FOREST = dict(lai=1.5, stand_density=.02, height=3., width=3., hspot=.025,
              optical_crown_geometry='clumping_dual', optical_crown_source='diffusion_profile')


class OpticalUrbanProfileTests(unittest.TestCase):
    def test_sunlit_wall_and_shaded_wall_have_different_brf(self):
        lit = urban_reflectance(CITY, 75., 0., 30., 0., .22, .3, .18, 0., 0.)
        shade = urban_reflectance(CITY, 75., 180., 30., 0., .22, .3, .18, 0., 0.)
        perpendicular = urban_reflectance(CITY, 75., 90., 30., 0., .22, .3, .18, 0., 0.)
        self.assertGreater(lit, shade+.05)
        self.assertGreater(lit, perpendicular+.05)

    def test_city_rotation_covariance(self):
        baseline = urban_reflectance(CITY, 60., 25., 30., 0., .22, .3, .18, 0., 0.)
        rotated = urban_reflectance(dict(CITY, urban_building_azimuth=73.),
                                    60., 98., 30., 73., .22, .3, .18, 0., 0.)
        self.assertAlmostEqual(baseline, rotated, places=12)

    def test_white_city_under_isotropic_sky(self):
        settings = dict(CITY, shortwave_diffuse_fraction=1.)
        for view, azimuth in ((0., 0.), (45., 90.), (75., 180.)):
            value = urban_reflectance(settings, view, azimuth, 30., 0., 1., 1., 1., 0., 0.)
            self.assertAlmostEqual(value, 1., places=10)

    def test_urban_reciprocity_and_escape_probability(self):
        for shape in ((10., 10., 20., .0015, 0.), (12., 4., 15., .002, 35.)):
            matrix, sky, walls = wall_exchange(shape, 0.)
            areas = np.r_[1-shape[0]*shape[1]*shape[3], walls]
            exchange = areas[:, None]*matrix
            np.testing.assert_allclose(exchange, exchange.T, atol=1e-10)
            np.testing.assert_allclose(matrix.sum(axis=1)+sky, 1., atol=1e-10)
            self.assertGreaterEqual(matrix.min(), 0.)

    def test_black_and_dark_city(self):
        self.assertEqual(urban_reflectance(CITY, 60., 0., 30., 0., 0., 0., 0., 0., 0.), 0.)
        settings = dict(CITY, shortwave_sunlit_surface=False)
        self.assertEqual(urban_reflectance(settings, 60., 0., 30., 0., .22, .3, .18, 0., 0.), 0.)

    def test_illumination_superposition(self):
        beam = urban_reflectance(CITY, 60., 0., 30., 0., .22, .3, .18, 0., 0.)
        sky = urban_reflectance(dict(CITY, shortwave_diffuse_fraction=1.), 60., 0., 30., 0., .22, .3, .18, 0., 0.)
        mixed = urban_reflectance(dict(CITY, shortwave_diffuse_fraction=.4), 60., 0., 30., 0., .22, .3, .18, 0., 0.)
        self.assertAlmostEqual(mixed, .6*beam+.4*sky, places=12)

    def test_no_buildings_recovers_forest_without_extra_cover_weight(self):
        settings = dict(CITY, **FOREST)
        settings.update(urban_density=0., urban_vegetation_fraction=1.)
        actual = urban_reflectance(settings, 45., 0., 30., 0., .22, .3, .18, .06, .03)
        expected = crown_reflectance(settings, 30., 45., 0., .06, .03, .18)
        self.assertEqual(actual, expected)

    def test_profile_conservative_limit_recovers_uniform_source(self):
        shape = (1.5, .02, 3., 3., .025)
        _, attenuation, upward, downward = source_profile(shape, .5, .5, 12)
        self.assertEqual(attenuation, 0.)
        self.assertAlmostEqual(upward, geometry(shape)[0], places=12)
        self.assertAlmostEqual(upward, downward, places=12)
        for view in (0., 30., 75.):
            new = crown_reflectance(FOREST, 30., view, 90., .5, .5, 1.)
            old = crown_reflectance(dict(FOREST, optical_crown_source='uniform'), 30., view, 90., .5, .5, 1.)
            self.assertAlmostEqual(new, old, places=12)

    def test_profile_energy_and_no_crown_spatial_integral(self):
        with patch('base.optical_crown.ellipsoid_grid_crown_FRT', side_effect=AssertionError('spatial crown grid')):
            for lai in (.1, 1.5, 6.):
                for solar in (0., 30., 75.):
                    for rho, tau, soil in ((.06, .03, .18), (.45, .42, .3), (.5, .5, 1.)):
                        _, budget = transport((lai, .02, 3., 3., .025), solar, rho, tau, soil,
                                               12, 'clumping_dual', 'diffusion_profile')
                        self.assertLess(abs(budget['balance_residual']), 1e-10)
                        if rho+tau == 1:
                            self.assertAlmostEqual(budget['reflected'], 1., places=10)

    def test_profile_does_not_enter_thermal_observations(self):
        config = dict(lai=1.5, hspot=.025, crown_density=.02, crown_rad_a=3., crown_rad_b=3.,
                      thermal_crown_geometry='clumping_dual', terrain_height=20., terrain_radius=30.,
                      terrain_density=.00015, building_density=.0015)
        spectra = dict(soil_e=.95, leaf_e=.985, terrain_e=.95, roof_e=.95, wall_e=.93, street_e=.955)
        with patch('base.optical_crown_profile.source_profile', side_effect=AssertionError('optical profile in thermal')):
            for surface in ('plane', 'slope', 'terrain', 'urban'):
                settings = dict(config, thermal_crown_geometry='legacy' if surface == 'urban' else 'clumping_dual')
                result = run._surface_result(settings, surface, 'crown', 10.5, (30., 0., 30., 0.), spectra)
                self.assertTrue(np.isfinite(result).all())


if __name__ == '__main__':
    unittest.main()
