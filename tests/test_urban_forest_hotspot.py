"""Urban forest two-scale geometry, physical limits and optical dispatch."""
import unittest
from unittest.mock import patch

import numpy as np

import run
from base.optical_urban import forest_parameters, urban_reflectance
from base.physicsF import planck
from base.thermal_crown import component_fractions
from base.urban_veg import Urban_Veg


CONFIG = dict(lai=1.5, hspot=.025, crown_density=.02, crown_rad_a=3., crown_rad_b=3.,
              building_length=10., building_width=10., building_height=20., building_density=.0015,
              thermal_crown_geometry='clumping_dual',
              street_sunlit_k=318.15, street_shaded_k=303.15,
              leaf_sunlit_k=304.65, leaf_shaded_k=304.65)
SPECTRA = dict(soil_e=.95, leaf_e=.985, terrain_e=.95, roof_e=.95, wall_e=.93, street_e=.955)
OPTICAL = dict(urban_building_length=10., urban_building_width=10., urban_building_height=20.,
               urban_density=.0015, urban_vegetation_fraction=1., lai=1.5,
               stand_density=.02, height=3., width=3., hspot=.025,
               optical_crown_geometry='clumping_dual', optical_crown_source='diffusion_profile')


def model(outer=1.):
    value = Urban_Veg()
    value.forestshape = np.array([1.5, .02, 3., 3., .025])
    value.set_vegetation_input('crown')
    value.set_hotspot_geometry('clumping_dual', .025, outer)
    return value


class UrbanForestHotspotTests(unittest.TestCase):
    def test_geometry_shared_with_thermal_adapter_and_coincident_limit(self):
        settings = dict(lai=1.5, stand_density=.02, height=3., width=3., hspot=.025,
                        thermal_crown_geometry='clumping_dual')
        view, azimuth = np.array([30., 25., 60.]), np.array([0., 90., 180.])
        actual = model().canopy_proportions(view, azimuth, 30., 0.)
        np.testing.assert_allclose(actual, component_fractions(settings, 30., view, azimuth), atol=1e-14)
        np.testing.assert_allclose(actual.sum(axis=0), 1., atol=1e-14)
        np.testing.assert_allclose(actual[[1, 3], 0], 0., atol=1e-14)
        self.assertGreaterEqual(actual.min(), 0.)

    def test_outer_hotspot_changes_thermal_radiance_and_preserves_marginals(self):
        first = model().canopy_proportions(25., 0., 30., 0.)
        off = model(0.).canopy_proportions(25., 0., 30., 0.)
        np.testing.assert_allclose(first[0]+first[1], off[0]+off[1], atol=1e-14)
        np.testing.assert_allclose(first[2]+first[3], off[2]+off[3], atol=1e-14)
        self.assertGreater(abs(first[0, 0]-off[0, 0]), 1e-4)
        active = run._surface_result(CONFIG, 'urban', 'crown', 10.5, (25., 0., 30., 0.), SPECTRA)
        disabled = run._surface_result(dict(CONFIG, thermal_hotspot_outer=0.),
                                       'urban', 'crown', 10.5, (25., 0., 30., 0.), SPECTRA)
        self.assertGreater(abs(active[0]-disabled[0]), 1e-5)

    def test_isothermal_has_no_hotspot_amplitude_and_blackbody_limit(self):
        config = dict(CONFIG)
        for material in ('soil', 'street', 'wall', 'roof', 'leaf'):
            for illumination in ('sunlit', 'shaded'):
                config[f'{material}_{illumination}_k'] = 300.
        for wave in (8., 10.5, 12.):
            for outer in (0., 1.):
                local = dict(config, thermal_hotspot_outer=outer)
                radiation, temperature = run._surface_result(
                    local, 'urban', 'crown', wave, (25., 0., 30., 0.), dict.fromkeys(SPECTRA, 1.))
                self.assertAlmostEqual(radiation, planck(wave, 300.), places=11)
                self.assertAlmostEqual(temperature, 300., places=10)
            a = run._surface_result(config, 'urban', 'crown', wave, (25., 0., 30., 0.), SPECTRA)
            b = run._surface_result(dict(config, thermal_hotspot_outer=0.),
                                   'urban', 'crown', wave, (25., 0., 30., 0.), SPECTRA)
            np.testing.assert_allclose(a, b, atol=1e-11, rtol=0)

    def test_empty_forest_and_no_building_limits(self):
        for empty in ({'lai': 0.}, {'crown_density': 0.}):
            active = run._surface_result(dict(CONFIG, **empty), 'urban', 'crown',
                                        10.5, (25., 0., 30., 0.), SPECTRA)
            legacy = run._surface_result(dict(CONFIG, **empty, thermal_crown_geometry='legacy'),
                                        'urban', 'crown', 10.5, (25., 0., 30., 0.), SPECTRA)
            np.testing.assert_allclose(active, legacy, atol=1e-11, rtol=0)
        config = dict(CONFIG, building_density=0., soil_sunlit_k=318.15, soil_shaded_k=303.15)
        urban = run._surface_result(config, 'urban', 'crown', 10.5, (25., 0., 30., 0.), SPECTRA)
        plane = run._surface_result(config, 'plane', 'crown', 10.5, (25., 0., 30., 0.),
                                    dict(SPECTRA, soil_e=SPECTRA['street_e']))
        np.testing.assert_allclose(urban, plane, atol=1e-11, rtol=0)

    def test_optical_outer_hotspot_is_already_active_and_not_duplicated(self):
        self.assertEqual(forest_parameters(OPTICAL)[3], (.025, 1.))
        for rho, tau, soil in ((.06, .03, .18), (.45, .42, .30)):
            original = urban_reflectance(OPTICAL, 25., 0., 30., 0., .22, .3, soil, rho, tau)
            explicit = urban_reflectance(dict(OPTICAL, optical_hotspot_inner=.025, optical_hotspot_outer=1.),
                                        25., 0., 30., 0., .22, .3, soil, rho, tau)
            off = urban_reflectance(dict(OPTICAL, optical_hotspot_outer=0.),
                                   25., 0., 30., 0., .22, .3, soil, rho, tau)
            self.assertEqual(original, explicit)
            self.assertGreater(abs(original-off), 1e-4)

    def test_parameterized_geometry_and_invalid_inputs(self):
        with patch('base.optical_crown.ellipsoid_grid_crown_FRT', side_effect=AssertionError('volume grid')):
            run._surface_result(CONFIG, 'urban', 'crown', 10.5, (24.781, 11.23, 30., 0.), SPECTRA)
        for mode, width in (('volume', 1.), ('typo', 1.), ('clumping_dual', -1.), ('clumping_dual', np.nan)):
            with self.assertRaises(ValueError):
                model().set_hotspot_geometry(mode, outer=width)


if __name__ == '__main__':
    unittest.main()
