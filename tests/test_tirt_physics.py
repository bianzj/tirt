"""Physical invariants and regression cases for TiRT surface/canopy coupling."""

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import angle_run
import run
from base.crown import Crown
from base.hom import Hom
from base.row import Row
from base.physicsF import planck, slope1


COMPONENTS = ("soil_e", "leaf_e", "terrain_e", "roof_e", "wall_e", "street_e")
SURFACES = ("plane", "slope", "terrain", "urban")
VEGETATIONS = ("bare", "hom", "row", "crown")
GEOMETRY = (40.0, 60.0, 35.0, 10.0)
CONFIG = {"lai": "1.5", "hspot": "0.15"}
BLACKBODY_RADIANCE = {8.0: 9.076384487115847, 10.5: 9.78997743993481, 12.0: 8.9600552403498}


def spectra(emissivity=0.95):
    return {name: emissivity for name in COMPONENTS}


def isothermal_config():
    result = CONFIG.copy()
    for material in ("soil", "leaf", "terrain", "roof", "wall", "street"):
        for lighting in ("sunlit", "shaded"):
            result[f"{material}_{lighting}_k"] = "300"
    return result


def direct_model(kind, lai=1.5, hotspot=0.15, row_gap=0.5, density=0.02):
    model = {"hom": Hom, "row": Row, "crown": Crown}[kind]()
    zenith = np.array([30.0])
    zero = np.array([0.0])
    if kind == "hom":
        model.set_structure(lai, hotspot)
        model.set_angle(zenith, zenith, zero)
    elif kind == "row":
        model.set_structure(lai, hotspot, 0.5, row_gap, 1.0)
        model.set_angle(zenith, zenith, zero, zero, zero)
    else:
        model.set_structure(lai, hotspot, density, 3.0, 3.0)
        model.set_angle(zenith, zenith, zero)
    model.set_thermal(300.0, 300.0, 300.0, 300.0)
    return model


class TiRTPhysicsTests(unittest.TestCase):
    def test_direct_model_multiband_blackbody(self):
        for kind in ("hom", "row", "crown"):
            model = direct_model(kind)
            for wavelength, expected_radiance in BLACKBODY_RADIANCE.items():
                with self.subTest(kind=kind, wavelength=wavelength):
                    model.set_optical(wavelength, 1.0, 1.0)
                    self.assertAlmostEqual(float(model.run(ifradiance=1)[0]), expected_radiance, places=9)
                    self.assertAlmostEqual(float(model.run()[0]), 300.0, places=9)

    def test_all_surface_canopy_blackbodies_at_all_bands(self):
        for runner in (run, angle_run):
            for surface in SURFACES:
                for vegetation in VEGETATIONS:
                    for wavelength, expected_radiance in BLACKBODY_RADIANCE.items():
                        with self.subTest(runner=runner.__name__, surface=surface, vegetation=vegetation, wavelength=wavelength):
                            radiance, temperature = runner._surface_result(
                                isothermal_config(), surface, vegetation, wavelength, GEOMETRY, spectra(1.0))
                            self.assertAlmostEqual(radiance, expected_radiance, places=9)
                            self.assertAlmostEqual(temperature, 300.0, places=9)

    def test_terrain_canopy_emissivity_scales_emission(self):
        for vegetation in ("hom", "row", "crown"):
            for emissivity in (0.1, 0.95):
                with self.subTest(vegetation=vegetation, emissivity=emissivity):
                    radiance, _ = run._surface_result(
                        isothermal_config(), "terrain", vegetation, 10.5, GEOMETRY, spectra(emissivity))
                    if vegetation == 'crown':
                        self.assertGreaterEqual(radiance, emissivity * BLACKBODY_RADIANCE[10.5])
                        self.assertLessEqual(radiance, BLACKBODY_RADIANCE[10.5])
                    else:
                        self.assertAlmostEqual(radiance, emissivity * BLACKBODY_RADIANCE[10.5], places=9)

    def test_bare_terrain_uses_soil_and_terrain_materials(self):
        baseline = run._surface_result(CONFIG, "terrain", "bare", 10.5, GEOMETRY, spectra())[0]
        roof = spectra() | {"roof_e": 0.3}
        self.assertAlmostEqual(run._surface_result(CONFIG, "terrain", "bare", 10.5, GEOMETRY, roof)[0], baseline, places=12)
        for component in ("soil_e", "terrain_e"):
            changed = spectra() | {component: 0.3}
            self.assertNotAlmostEqual(run._surface_result(CONFIG, "terrain", "bare", 10.5, GEOMETRY, changed)[0], baseline, places=6)

    def test_canopy_models_are_really_dispatched(self):
        for surface in SURFACES:
            values = [run._surface_result(CONFIG, surface, kind, 10.5, GEOMETRY, spectra())[0]
                      for kind in ("hom", "row", "crown")]
            with self.subTest(surface=surface):
                self.assertEqual(len({round(value, 6) for value in values}), 3)

    def test_row_dimensions_and_bearing_affect_every_surface(self):
        for surface in SURFACES:
            for key, alternatives in (("row_width", (0.25, 1.0)), ("row_height", (1.0, 2.0)), ("row_azimuth", (0.0, 90.0))):
                values = [run._surface_result(CONFIG | {key: str(value)}, surface, "row", 10.5, GEOMETRY, spectra())[0]
                          for value in alternatives]
                with self.subTest(surface=surface, key=key):
                    self.assertGreater(abs(values[0] - values[1]), 1.0e-5)

    def test_zero_slope_matches_corresponding_plane(self):
        for kind in VEGETATIONS:
            config = CONFIG | {"slope_angle": "0", "slope_aspect": "45"}
            plane = run._surface_result(config, "plane", kind, 10.5, GEOMETRY, spectra())[0]
            slope = run._surface_result(config, "slope", kind, 10.5, GEOMETRY, spectra())[0]
            with self.subTest(kind=kind):
                self.assertAlmostEqual(plane, slope, places=6)

    def test_zero_lai_direct_api_is_sunlit_bare_soil(self):
        for kind in ("hom", "row", "crown"):
            model = direct_model(kind, lai=0.0)
            model.set_optical(10.5, 0.95, 0.985)
            with self.subTest(kind=kind):
                self.assertAlmostEqual(float(model.run(ifradiance=1)[0]), 0.95 * BLACKBODY_RADIANCE[10.5], places=9)

    def test_zero_crown_density_direct_api_has_no_leaf_scatter(self):
        model = direct_model("crown", density=0.0)
        model.set_optical(10.5, 0.95, 0.985)
        self.assertAlmostEqual(float(model.run(ifradiance=1)[0]), 0.95 * BLACKBODY_RADIANCE[10.5], places=9)

    def test_zero_hotspot_direct_apis_remain_finite(self):
        for kind in ("hom", "row", "crown"):
            model = direct_model(kind, hotspot=0.0)
            model.set_optical(10.5, 0.95, 0.985)
            with self.subTest(kind=kind):
                self.assertTrue(np.all(np.isfinite(model.run())))

    def test_zero_row_gap_matches_homogeneous_canopy(self):
        hom, row = direct_model("hom"), direct_model("row", row_gap=0.0)
        for model in (hom, row):
            model.set_optical(10.5, 0.95, 0.985)
        np.testing.assert_allclose(hom.run(ifradiance=1), row.run(ifradiance=1), atol=1.0e-12)

    def test_backlit_slope_zero_lai_uses_shaded_soil(self):
        config = CONFIG | {"lai": "0", "slope_angle": "60", "slope_aspect": "0"}
        geometry = (0.0, 0.0, 60.0, 180.0)
        expected = 0.95 * float(planck(10.5, 303.15))
        for kind in VEGETATIONS:
            with self.subTest(kind=kind):
                radiance, _ = run._surface_result(config, "slope", kind, 10.5, geometry, spectra())
                self.assertAlmostEqual(radiance, expected, places=9)

    def test_empty_buildings_and_terrain_reduce_to_plane(self):
        for kind in VEGETATIONS:
            plane = run._surface_result(CONFIG, "plane", kind, 10.5, GEOMETRY, spectra())
            for surface, key in (("urban", "building_density"), ("terrain", "terrain_density")):
                with self.subTest(kind=kind, surface=surface):
                    result = run._surface_result(CONFIG | {key: "0"}, surface, kind, 10.5, GEOMETRY, spectra())
                    np.testing.assert_allclose(result, plane, atol=1.0e-12)

    def test_nadir_azimuth_is_invariant(self):
        for surface in SURFACES:
            for kind in VEGETATIONS:
                values = [run._surface_result(CONFIG, surface, kind, 10.5, (0.0, azimuth, 30.0, 0.0), spectra())[0]
                          for azimuth in (0.0, 90.0, 180.0, 270.0)]
                with self.subTest(surface=surface, kind=kind):
                    np.testing.assert_allclose(values, [values[0]] * 4, atol=1.0e-9)

    def test_horizon_blackbody_limit_and_explicit_invalid_kernel(self):
        for surface in SURFACES:
            for kind in VEGETATIONS:
                with self.subTest(surface=surface, kind=kind):
                    result = run._surface_result(isothermal_config(), surface, kind, 10.5, (90.0, 0.0, 30.0, 0.0), spectra(1.0))
                    if surface == "slope":
                        # A horizontal plane is tangent to this view and is
                        # excluded by the explicit slope visibility rule.
                        self.assertTrue(np.isnan(result[1]))
                    else:
                        self.assertAlmostEqual(result[1], 300.0, places=8)
        # The legacy bare-terrain scattering quadrature becomes nonphysical
        # near the horizon. It must fail explicitly, never switch to a plane.
        with self.assertRaisesRegex(FloatingPointError, "terrain\\+bare"):
            run._surface_result(CONFIG, "terrain", "bare", 10.5, (89.9, 0.0, 30.0, 0.0), spectra())

    def test_slope_angle_transform_has_exact_nadir_and_azimuth(self):
        zenith, azimuth = slope1(np.array([0.0, 45.0, 45.0, 120.0]),
                                 np.array([270.0, 45.0, 180.0, 0.0]), 0.0, 0.0)
        np.testing.assert_allclose(zenith, [0.0, 45.0, 45.0, 60.0], atol=1.0e-5)
        np.testing.assert_allclose(azimuth, [0.0, 45.0, 180.0, 0.0], atol=1.0e-12)

    def test_numeric_failure_is_not_replaced_by_bare_soil(self):
        with patch.object(Hom, "run", return_value=np.array([np.nan])):
            with self.assertRaisesRegex(FloatingPointError, "plane\\+hom"):
                run._vegetation_result(CONFIG, "hom", 10.5, GEOMETRY, spectra())


if __name__ == "__main__":
    unittest.main()
