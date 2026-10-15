"""Forest terrain dispatch, energy limits and angular transformations."""
import math
import unittest
from unittest.mock import patch

import numpy as np

import run
from base.optical import reflectance
from base.optical_crown import crown_reflectance
from base.optical_terrain import terrain_crown_reflectance, exchange_geometry
from base.terrain_rt import _mountain_visibility
from base.terrain_forest import TerrainForest


SETTINGS = dict(surface_model='terrain', terrain_type=2, vegetation_model='crown',
                surface_type=1, distribution=0, lai=1.5, hspot=.025, height=3., width=3.,
                stand_density=.02, terrain_height=20., terrain_radius=30., terrain_density=.00015,
                optical_crown_geometry='clumping_dual', optical_crown_source='uniform',
                optical_crown_quadrature=4, wavelengths=[670, 860],
                leaf_reflectance_spectrum=[.06, .45], leaf_transmittance_spectrum=[.03, .42],
                soil_reflectance_spectrum=[.18, .3])
CONFIG = dict(lai=1.5, hspot=.025, crown_density=.02, crown_rad_a=3., crown_rad_b=3.,
              crown_height=3., thermal_crown_geometry='clumping_dual',
              terrain_height=20., terrain_radius=30., terrain_density=.00015,
              soil_sunlit_k=318.15, soil_shaded_k=303.15,
              leaf_sunlit_k=306.15, leaf_shaded_k=303.15)
SPECTRA = dict(soil_e=.95, leaf_e=.985, terrain_e=.95)


class ForestTerrainTests(unittest.TestCase):
    def test_cone_normal_angle_matches_actual_geometry(self):
        _, _, details = _mountain_visibility(((20., 30., .00015),), 30., 30., 0.)
        self.assertAlmostEqual(details[0][0], math.atan2(20., 30.), places=14)

    def test_no_terrain_matches_plane(self):
        settings = dict(SETTINGS, terrain_density=0.)
        actual = terrain_crown_reflectance(settings, 30., 0., 25., 0., .45, .42, .3)
        expected = crown_reflectance(dict(settings, surface_model='plane'), 30., 25., 0., .45, .42, .3)
        self.assertEqual(actual, expected)
        actual = run._surface_result(dict(CONFIG, terrain_density=0.), 'terrain', 'crown', 10.5,
                                     (25., 0., 30., 0.), SPECTRA)
        expected = run._surface_result(CONFIG, 'plane', 'crown', 10.5, (25., 0., 30., 0.), SPECTRA)
        np.testing.assert_array_equal(actual, expected)

    def test_terrain_optics_uses_dual_hotspots_without_spatial_integral(self):
        with patch('base.radiation._component_proportions', side_effect=AssertionError('legacy dispatch')), \
             patch('base.optical_crown.ellipsoid_grid_crown_FRT', side_effect=AssertionError('spatial crown grid')):
            baseline = reflectance(SETTINGS, 860., 25., 0., 30., 0.)
            changed = reflectance(dict(SETTINGS, optical_hotspot_outer=.02), 860., 25., 0., 30., 0.)
            self.assertGreater(abs(baseline-changed), 1e-4)

    def test_terrain_exchange_and_illumination_are_conservative(self):
        matrix = exchange_geometry(((20., 30., .00015),))
        self.assertTrue((matrix >= 0).all())
        self.assertTrue((matrix.sum(axis=1) <= 1).all())
        beam = terrain_crown_reflectance(SETTINGS, 30., 0., 25., 0., .45, .42, .3)
        sky = terrain_crown_reflectance(SETTINGS, 30., 0., 25., 0., .45, .42, .3, 1.)
        mixed = terrain_crown_reflectance(SETTINGS, 30., 0., 25., 0., .45, .42, .3, .3)
        self.assertAlmostEqual(mixed, .7*beam+.3*sky, places=12)
        shaded = terrain_crown_reflectance(SETTINGS, 30., 0., 25., 0., .45, .42, .3, 0., False)
        self.assertEqual(shaded, 0.)
        black = terrain_crown_reflectance(SETTINGS, 30., 0., 25., 0., 0., 0., 0.)
        self.assertEqual(black, 0.)
        with self.assertRaises(ValueError):
            terrain_crown_reflectance(dict(SETTINGS, terrain_density=.01), 30., 0., 25., 0., .45, .42, .3)

    def test_zero_slope_thermal_dual_matches_plane(self):
        geometry = (25., 0., 30., 0.)
        slope = run._surface_result(dict(CONFIG, slope_angle=0., slope_aspect=90.),
                                   'slope', 'crown', 10.5, geometry, SPECTRA)
        plane = run._surface_result(CONFIG, 'plane', 'crown', 10.5, geometry, SPECTRA)
        np.testing.assert_allclose(slope, plane, atol=1e-12)

    def test_thermal_outer_hotspot_is_active_on_slope_and_terrain(self):
        for surface in ('slope', 'terrain'):
            config = dict(CONFIG, slope_angle=30., slope_aspect=90.)
            baseline = run._surface_result(config, surface, 'crown', 10.5, (25., 0., 30., 0.), SPECTRA)[0]
            changed = run._surface_result(dict(config, thermal_hotspot_outer=.02), surface,
                                         'crown', 10.5, (25., 0., 30., 0.), SPECTRA)[0]
            self.assertGreater(abs(baseline-changed), 1e-4)

    def test_thermal_scatter_restored_and_blackbody_limit(self):
        model = TerrainForest()
        model.forestshape = np.array([1.5, .02, 3., 3., .025, 3.])
        model.set_hotspot_geometry('clumping_dual')
        model.set_structural_input([[20., 30., .00015]])
        model.set_angular_input(np.array([30.]), np.array([0.]), 30., 0.)
        model.set_spectral_input(.95, .985)
        direct = np.asarray(model.calculate_component_direct_emissivity(0, 2))
        total = model.calculate_effective_component_emissivity(2)[0]
        self.assertAlmostEqual(model.calculate_effective_component_emissivity(0)[0], total.sum())
        self.assertGreater(total.sum(), direct.sum())
        self.assertLessEqual(total.sum(), 1.)
        model.set_spectral_input(1., 1.)
        self.assertAlmostEqual(model.calculate_effective_component_emissivity(2).sum(), 1., places=12)


if __name__ == '__main__':
    unittest.main()
