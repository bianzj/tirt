"""Physical limits and entry-point checks for the thermal two-scale adapter."""
import unittest
from unittest.mock import patch

import numpy as np

import run
import angle_run
from base.brightness import vegetation_brightness
from base.crown import Crown
from base.io import read_defined
from base.physicsF import planck
from base.thermal_crown import component_fractions, thermal_emission


SETTINGS = dict(lai=1., stand_density=1/64, height=1.5, width=2., hspot=.2/3,
                soil_emissivity=.95, leaf_emissivity=.98, vegetation_model="crown",
                surface_type=1, distribution=1, leaf_reflectance=.015,
                leaf_transmittance=.005, soil_reflectance=.05)
TEMPERATURES = [315., 295., 305., 298.]
STATE = dict(zip(("Tss", "Tsh", "Tls", "Tlh"), TEMPERATURES), sza=30., saa=0.)
MODES = ("clumping_dual", "volume")


class ThermalCrownTests(unittest.TestCase):
    def test_probabilities_marginals_and_coincident_limit(self):
        for mode in MODES:
            s = dict(SETTINGS, thermal_crown_geometry=mode)
            for solar in (0., 30., 60.):
                for view in (0., 30., 60., 85.):
                    p = component_fractions(s, solar, view, [0., 90., 180.])
                    self.assertTrue(np.isfinite(p).all())
                    self.assertGreaterEqual(p.min(), -1e-14)
                    np.testing.assert_allclose(p.sum(axis=0), 1., atol=1e-13)
                    np.testing.assert_allclose(p[0]+p[1], (p[0]+p[1])[0], atol=1e-13)
                p = component_fractions(s, solar, solar, 0.)
                np.testing.assert_allclose(p[[1, 3]], 0., atol=1e-13)

    def test_blackbody_isothermal_multiband(self):
        for mode in MODES:
            s = dict(SETTINGS, thermal_crown_geometry=mode, soil_emissivity=1., leaf_emissivity=1.)
            for wave in (3.7, 8., 10.5, 12.):
                r, t, _ = thermal_emission(s, [300.]*4, wave, 30., [0., 30., 60.], [0., 0., 180.])
                np.testing.assert_allclose(r, planck(wave, 300.), rtol=1e-13)
                np.testing.assert_allclose(t, 300., atol=1e-10)

    def test_no_sun_shade_contrast_has_no_azimuth_hotspot(self):
        for mode in MODES:
            s = dict(SETTINGS, thermal_crown_geometry=mode)
            for temps in ([300.]*4, [310., 310., 300., 300.]):
                r, _, _ = thermal_emission(s, temps, 10.5, 30., 30., [0., 90., 180.])
                np.testing.assert_allclose(r, r[0], atol=1e-13)

    def test_sparse_reflective_crowns_do_not_exceed_blackbody_budget(self):
        for mode in MODES:
            for density in (.0001, 1/64, .1):
                s = dict(SETTINGS, thermal_crown_geometry=mode, lai=6., stand_density=density,
                         soil_emissivity=1., leaf_emissivity=.5)
                r, t, _ = thermal_emission(s, [300.]*4, 10.5, 30., 30., [0., 90., 180.])
                self.assertLessEqual(r.max(), planck(10.5, 300.)+1e-12)
                self.assertLessEqual(t.max(), 300.+1e-10)
                np.testing.assert_allclose(r, r[0], atol=1e-13)

    def test_hotspot_sign_follows_temperature_contrast(self):
        for mode in MODES:
            s = dict(SETTINGS, thermal_crown_geometry=mode)
            hot = thermal_emission(s, TEMPERATURES, 10.5, 30., 30., [0., 90.])[1]
            cold = thermal_emission(s, [295., 315., 298., 305.], 10.5, 30., 30., [0., 90.])[1]
            self.assertGreater(hot[0], hot[1])
            self.assertLess(cold[0], cold[1])

    def test_night_ignores_sunlit_temperatures(self):
        for mode in MODES:
            s = dict(SETTINGS, thermal_crown_geometry=mode)
            a = thermal_emission(s, TEMPERATURES, 10.5, 100., 30., [0., 90.], sunlit=False)
            b = thermal_emission(s, [400., 295., 400., 298.], 10.5, 100., 30., [0., 90.], sunlit=False)
            np.testing.assert_array_equal(a[0], b[0])
            np.testing.assert_array_equal(a[2][[0, 2]], 0.)
            np.testing.assert_allclose(a[0], a[0][0], atol=1e-13)

    def test_bare_and_zero_density_limits(self):
        for mode in MODES:
            for empty in ({"lai": 0.}, {"stand_density": 0.}):
                s = dict(SETTINGS, thermal_crown_geometry=mode, **empty)
                for sunlit, temp in ((True, 315.), (False, 295.)):
                    r, _, _ = thermal_emission(s, TEMPERATURES, 10.5, 30., 60., 90., sunlit)
                    np.testing.assert_allclose(r, .95*planck(10.5, temp), rtol=1e-13)

    def test_widths_are_independent_and_optical_settings_do_not_leak(self):
        s = dict(SETTINGS, thermal_crown_geometry="clumping_dual")
        p = component_fractions(s, 30., 25., 0.)
        for key in ("thermal_hotspot_inner", "thermal_hotspot_outer"):
            off = component_fractions(dict(s, **{key: 0.}), 30., 25., 0.)
            self.assertGreater(np.max(np.abs(off-p)), 1e-5)
            np.testing.assert_allclose(off[0]+off[1], p[0]+p[1], atol=1e-13)
        np.testing.assert_array_equal(component_fractions(dict(s, optical_hotspot_inner=0., optical_hotspot_outer=0.), 30., 25., 0.), p)

    def test_parameterized_path_has_no_volume_quadrature(self):
        with patch("base.optical_crown.proportion_bidirectional_crown_voxel_one", side_effect=AssertionError("volume called")):
            component_fractions(dict(SETTINGS, thermal_crown_geometry="clumping_dual"), 31.123, 25.123, 4.123)

    def test_all_plane_entry_points_match(self):
        for mode in MODES:
            s = dict(SETTINGS, thermal_crown_geometry=mode)
            expected = thermal_emission(s, TEMPERATURES, 10.5, 30., 25., 0.)[0][0]
            self.assertAlmostEqual(vegetation_brightness(s, STATE, 10.5, 25., 0.)[0], expected, places=12)
            model = Crown()
            model.set_structure(1., .2/3, 1/64, 1.5, 2.)
            model.set_hotspot_geometry(mode)
            model.set_angle(np.array([25.]), np.array([30.]), np.array([0.]))
            model.set_optical(10.5, .95, .98)
            model.set_thermal(*TEMPERATURES)
            self.assertAlmostEqual(model.run(1)[0], expected, places=12)
            config = dict(lai=1., hspot=.2/3, crown_density=1/64, crown_rad_b=1.5, crown_rad_a=2., thermal_crown_geometry=mode)
            config.update(zip(("soil_sunlit_k", "soil_shaded_k", "leaf_sunlit_k", "leaf_shaded_k"), TEMPERATURES))
            for entry in (run, angle_run):
                actual = entry._surface_result(config, "plane", "crown", 10.5, (25., 0., 30., 0.), dict(soil_e=.95, leaf_e=.98))[0]
                self.assertAlmostEqual(actual, expected, places=12)

    def test_settings_reader_propagates_mode(self):
        s = read_defined(overrides={"thermal_crown_geometry": "clumping_dual", "thermal_hotspot_inner": "0.08", "thermal_hotspot_outer": "1.2"})
        self.assertEqual(s["thermal_crown_geometry"], "clumping_dual")
        self.assertEqual(float(s["thermal_hotspot_inner"]), .08)
        self.assertEqual(float(s["thermal_hotspot_outer"]), 1.2)

    def test_invalid_inputs_fail_explicitly(self):
        for change in ({"thermal_crown_geometry": "typo"}, {"thermal_hotspot_outer": -1.},
                       {"thermal_crown_geometry": "volume", "thermal_hotspot_inner": .1}):
            s = dict(SETTINGS, thermal_crown_geometry="clumping_dual")
            s.update(change)
            with self.assertRaises(ValueError):
                component_fractions(s, 30., 25., 0.)


if __name__ == "__main__":
    unittest.main()
