"""Terrain-energy controls use the real plane radiation implementation."""
import importlib
import math
import unittest
import numpy as np
from base import tirteb_radiation as radiation
from base.solar import solar_angles

try:
    terrain = importlib.import_module('base.terrain_energy')
except ModuleNotFoundError:
    terrain = None

PLANE_SW = getattr(radiation, '_plane_prepare_shortwave', radiation.prepare_shortwave)
PLANE_LW = getattr(radiation, '_plane_prepare_longwave', radiation.prepare_longwave)
PLANE_WEIGHTS = getattr(radiation, '_plane_flux_weights', radiation.flux_weights)


class TerrainEnergyTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(terrain, 'terrain must participate in the energy budget')
        self.settings = dict(surface_model='terrain', terrain_type=1, terrain_slope=40.,
                             terrain_aspect=0., vegetation_model='bare', lai=0.,
                             soil_reflectance=0., leaf_reflectance=0.1,
                             leaf_transmittance=0.1, soil_emissivity=1., leaf_emissivity=1.,
                             shortwave_diffuse_fraction=0.2,
                             shortwave_enforce_solar_limit=0)
        self.meteo = dict(shortwave=962., longwave=350.)
        self.state = dict(Tss=300., Tsh=300., Tls=300., Tlh=300.)

    def sw(self, settings=None, sza=35., saa=160.):
        return terrain.prepare_shortwave(settings or self.settings, self.meteo, sza, saa,
                                         plane_shortwave=PLANE_SW, plane_weights=PLANE_WEIGHTS)

    def weights(self, settings=None, sza=35., saa=160.):
        return np.asarray(terrain.flux_weights(settings or self.settings, sza, saa,
                                              plane_weights=PLANE_WEIGHTS))

    def test_audit_north_slope_uses_horizontal_direct_projection(self):
        sza, saa = solar_angles(2024, 214, 214.5556-8/24, 39.9, 116.4)
        self.assertAlmostEqual(sza, 25.314277165002792)
        result = self.sw(sza=sza, saa=saa)
        self.assertAlmostEqual(result['shortwave_direct_incident']/769.6, 0.5157735712736135)

    def test_backfacing_slope_retains_isotropic_diffuse(self):
        result = self.sw(sza=80., saa=180.)
        self.assertAlmostEqual(result['shortwave_direct_incident'], 0.)
        self.assertAlmostEqual(result['diffuseVrad_soil_shaded'],
                               962.*0.2*(1+math.cos(math.radians(40)))/2)
        self.assertEqual(self.weights(sza=80., saa=180.)[0], 0)

    def test_flat_slope_reduces_to_original_plane_budget(self):
        settings = dict(self.settings, terrain_slope=0., soil_reflectance=0.2)
        actual = self.sw(settings)
        expected = PLANE_SW(dict(settings, surface_model='plane'), self.meteo, 35., 160.)
        for key in ('directVrad_soil', 'diffuseVrad_soil', 'shortwave_direct_incident'):
            self.assertAlmostEqual(actual[key], expected[key])
        np.testing.assert_allclose(self.weights(settings), [1, 0, 0, 0])

    def test_zero_density_compound_reduces_to_plane(self):
        settings = dict(self.settings, terrain_type=2, terrain_density=0., terrain_height=20., terrain_radius=10.)
        actual = self.sw(settings)
        self.assertAlmostEqual(actual['shortwave_direct_incident'], 769.6)
        self.assertAlmostEqual(actual['diffuseVrad_soil_sunlit'], 192.4)
        np.testing.assert_allclose(self.weights(settings), [1, 0, 0, 0])

    def test_flat_absent_shade_keeps_plane_diffuse_and_longwave_limit(self):
        settings = dict(self.settings, terrain_slope=0.)
        actual = self.sw(settings)
        self.assertAlmostEqual(actual['diffuseVrad_soil_shaded'], 192.4)
        lw = terrain.prepare_longwave(settings, self.meteo, self.state, 35., 160.,
                                     plane_longwave=PLANE_LW, plane_weights=PLANE_WEIGHTS)
        self.assertAlmostEqual(lw['LW_soil_shaded'], 350.)

    def test_compound_canopy_preserves_lai_and_shortwave_energy_bounds(self):
        for kind in ('hom', 'row', 'crown'):
            settings = dict(self.settings, terrain_type=2, terrain_density=0.001,
                            terrain_height=20., terrain_radius=10., vegetation_model=kind,
                            lai=2., stand_density=0.02, height=3., width=3., soil_reflectance=0.25)
            result = self.sw(settings, sza=70.)
            weights = self.weights(settings, sza=70.)
            absorbed = 0.
            for index, name in ((0, 'soil'), (2, 'leaf')):
                absorbed += weights[index]*(result['directVrad_'+name]+result['diffuseVrad_'+name+'_sunlit'])
                absorbed += weights[index+1]*result['diffuseVrad_'+name+'_shaded']
            self.assertAlmostEqual(weights[2]+weights[3], 2.)
            self.assertGreater(absorbed, 0.)
            self.assertLessEqual(absorbed, 962.+1e-7)
            self.assertAlmostEqual(absorbed, 962.*(1-result['shortwave_albedo']))

    def test_compound_first_hits_conserve_horizontal_incident_power(self):
        settings = dict(self.settings, terrain_type=2, terrain_density=0.001,
                        terrain_height=20., terrain_radius=10.)
        for sza in (0., 35., 80.):
            result = self.sw(settings, sza=sza)
            weights = self.weights(settings, sza=sza)
            absorbed = weights[0]*(result['directVrad_soil']+result['diffuseVrad_soil_sunlit']) \
                       + weights[1]*result['diffuseVrad_soil_shaded']
            self.assertAlmostEqual(absorbed, 962., places=7)
            self.assertAlmostEqual(result['terrain_direct_intercepted_ground'], 769.6)
            self.assertGreater(weights[0]+weights[1], 1.)

    def test_compound_view_factors_are_reciprocal_and_sky_area_is_horizontal(self):
        settings = dict(self.settings, terrain_type=2, terrain_density=0.001,
                        terrain_height=20., terrain_radius=10.)
        layout = terrain.layout(settings, 45., 160.)
        areas, factors, sky = layout['area'], layout['factors'], layout['sky']
        np.testing.assert_allclose(areas[:, None]*factors, areas[None, :]*factors.T, atol=1e-12)
        np.testing.assert_allclose(factors.sum(axis=1)+sky, 1.)
        self.assertAlmostEqual(float(areas @ sky), 1.)
        self.assertTrue(np.all(factors >= 0) and np.all(sky >= 0))

    def test_equal_temperature_enclosure_has_zero_component_longwave_net(self):
        for terrain_type, lai in ((1, 0.), (2, 0.), (2, 2.)):
            settings = dict(self.settings, terrain_type=terrain_type, terrain_density=0.001,
                            terrain_height=20., terrain_radius=10., vegetation_model='hom', lai=lai,
                            soil_emissivity=0.95, leaf_emissivity=0.97)
            meteo = dict(self.meteo, longwave=5.67e-8*300**4)
            result = terrain.prepare_longwave(settings, meteo, self.state, 35., 160.,
                                              plane_longwave=PLANE_LW, plane_weights=PLANE_WEIGHTS)
            weights = self.weights(settings)
            for index, name in enumerate(('soil_sunlit','soil_shaded','leaf_sunlit','leaf_shaded')):
                if weights[index] > 1e-10:
                    expected = (0.95 if index<2 else 2*0.97)*5.67e-8*300**4
                    self.assertAlmostEqual(result['LW_'+name], expected, places=7)


if __name__ == '__main__':
    unittest.main()

