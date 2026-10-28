"""Independent leaf physiology and broadband/observation separation regressions."""
import math
from pathlib import Path
import tempfile
import unittest

from base.biochemical import suresist
from base.eb import air_density
from base.io import read_defined


ROOT = Path(__file__).resolve().parents[1]
METEO = {"pressure": 1013.25, "vapor_pressure": 15.0, "air_temperature_c": 25.0}
FIXTURE = tempfile.TemporaryDirectory()
COEFFICIENTS = Path(FIXTURE.name) / "empty_optical_coefficients.csv"
COEFFICIENTS.write_text("", encoding="utf-8")


def settings(**overrides):
    # Exercise scalar fallback optics; FLUSPECT's leaf model is outside these tests.
    return read_defined(overrides={"_input_base": str(ROOT),
                                   "optical_coefficients_file": str(COEFFICIENTS), **overrides})


def response(config, **changes):
    state = {"Tls": 298.15, "Tlh": 298.15, "cs": 380.0, "es": 15.0}
    state.update(changes)
    return suresist(config, METEO, state, {
        "directPrad_leaf": 900.0, "diffusePrad_leaf": 150.0,
        "diffusePrad_leaf_sunlit": 180.0, "diffusePrad_leaf_shaded": 120.0,
    })


class LeafPhysiologyTests(unittest.TestCase):
    def test_shaded_leaf_does_not_receive_direct_par(self):
        for method in ("farquhar", "ballberry"):
            with self.subTest(method=method):
                config = settings(stomatal_method=method)
                result = response(config)
                self.assertLess(result["rssunlit"], result["rsshaded"])
                state = {"Tls": 298.15, "Tlh": 298.15, "cs": 380.0, "es": 15.0}
                dark_shade = suresist(config, METEO, state, {
                    "directPrad_leaf": 900.0, "diffusePrad_leaf": 150.0,
                    "diffusePrad_leaf_sunlit": 180.0, "diffusePrad_leaf_shaded": 0.0,
                })
                self.assertGreater(dark_shade["rsshaded"], result["rsshaded"])
                self.assertAlmostEqual(dark_shade["rssunlit"], result["rssunlit"])

    def test_farquhar_dark_leaf_keeps_minimum_conductance(self):
        config = settings(stomatal_method="farquhar", stomatal_g0=0.02)
        state = {"Tls": 298.15, "Tlh": 298.15, "cs": 380.0, "es": 15.0}
        dark = suresist(config, METEO, state, {"directPrad_leaf": 0.0, "diffusePrad_leaf_sunlit": 0.0,
                                               "diffusePrad_leaf_shaded": 0.0})
        g0_resistance = air_density(METEO) / 28.96e-3 / 0.02
        self.assertAlmostEqual(dark["rssunlit_leaf"], g0_resistance, delta=1.0e-6 * g0_resistance)
        self.assertAlmostEqual(dark["rsshaded_leaf"], g0_resistance, delta=1.0e-6 * g0_resistance)
        lit = response(config)
        self.assertLess(lit["rssunlit_leaf"], g0_resistance)

    def test_each_temperature_affects_only_its_leaf(self):
        for method in ("farquhar", "ballberry"):
            with self.subTest(method=method):
                config = settings(stomatal_method=method)
                reference = response(config)
                hot_sun = response(config, Tls=313.15)
                hot_shade = response(config, Tlh=313.15)
                self.assertNotAlmostEqual(reference["rssunlit"], hot_sun["rssunlit"])
                self.assertAlmostEqual(reference["rsshaded"], hot_sun["rsshaded"])
                self.assertNotAlmostEqual(reference["rsshaded"], hot_shade["rsshaded"])
                self.assertAlmostEqual(reference["rssunlit"], hot_shade["rssunlit"])

    def test_resistance_is_converted_to_ground_area_bulk_units(self):
        for method in ("farquhar", "ballberry"):
            with self.subTest(method=method):
                one = response(settings(lai=1.0, stomatal_method=method))
                three = response(settings(lai=3.0, stomatal_method=method))
                self.assertAlmostEqual(one["rssunlit"], 3.0 * three["rssunlit"])
                self.assertAlmostEqual(one["rsshaded"], 3.0 * three["rsshaded"])

    def test_default_soil_resistance_uses_supplied_saturation(self):
        dry = response(settings(soil_moisture=0.0, satwater=0.6))["rss"]
        wet = response(settings(soil_moisture=0.6, satwater=0.6))["rss"]
        self.assertAlmostEqual(dry, math.exp(8.206))
        self.assertAlmostEqual(wet, math.exp(3.951))
        self.assertGreater(dry, wet)

    def test_explicit_constant_soil_resistance_is_preserved(self):
        for overrides in ({"soil_rss": 123.0}, {"rss": 123.0},
                          {"soil_rss": 123.0, "soil_rss_method": "constant"}):
            for water in (0.0, 0.6):
                self.assertEqual(response(settings(soil_moisture=water, satwater=0.6, **overrides))["rss"], 123.0)

    def test_explicit_sellers_mode_overrides_constant_value(self):
        result = response(settings(soil_moisture=0.3, satwater=0.6, soil_rss=123.0, soil_rss_method="sellers"))
        self.assertAlmostEqual(result["rss"], math.exp(6.0785))


class BroadbandTests(unittest.TestCase):
    def test_observation_sampling_does_not_change_energy_parameters(self):
        configs = [settings(wavelengths="8;10.5", leaf_reflectance="0.1;0.02", soil_reflectance="0.15;0.03",
                            roof_reflectance="0.1;0.02", wall_reflectance="0.15;0.03", street_reflectance="0.2;0.05"),
                   settings(wavelengths="8;8.1;8.2;10.5", leaf_reflectance="0.1;0.1;0.1;0.02", soil_reflectance="0.15;0.15;0.15;0.03",
                            roof_reflectance="0.1;0.1;0.1;0.02", wall_reflectance="0.15;0.15;0.15;0.03", street_reflectance="0.2;0.2;0.2;0.05"),
                   settings(wavelengths="0.67")]
        for key in ("leaf_emissivity", "leaf_reflectance", "leaf_transmittance", "soil_emissivity",
                    "urban_emissivity_roof", "urban_emissivity_wall", "urban_emissivity_street"):
            self.assertEqual(configs[0][key], configs[1][key], key)
            self.assertEqual(configs[0][key], configs[2][key], key)
        self.assertEqual(configs[0]["observation_bands"][0]["leaf_reflectance"], 0.1)
        self.assertAlmostEqual(configs[0]["observation_bands"][1]["soil_reflectance"], 0.03)

    def test_explicit_broadband_emissivities_do_not_replace_observation_optics(self):
        config = settings(wavelengths="10.5", leaf_reflectance="0.07", soil_reflectance="0.12",
                          broadband_leaf_emissivity=0.91, broadband_soil_emissivity=0.88,
                          broadband_urban_emissivity_roof=0.81, broadband_urban_emissivity_wall=0.82,
                          broadband_urban_emissivity_street=0.83)
        self.assertEqual(config["leaf_emissivity"], 0.91)
        self.assertEqual(config["soil_emissivity"], 0.88)
        self.assertEqual(config["urban_emissivity_roof"], 0.81)
        self.assertEqual(config["urban_emissivity_wall"], 0.82)
        self.assertEqual(config["urban_emissivity_street"], 0.83)
        self.assertEqual(config["observation_bands"][0]["leaf_reflectance"], 0.07)
        self.assertEqual(config["observation_bands"][0]["soil_reflectance"], 0.12)

    def test_nonphysical_broadband_emissivity_is_rejected(self):
        for value in (-0.1, 1.01, float("nan")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                settings(broadband_leaf_emissivity=value)

    def test_supplied_timestep_and_saturation_are_exposed(self):
        config = settings(dtime=900.0, satwater=0.6, soil_moisture=0.5)
        self.assertEqual(config["dtime"], 900.0)
        self.assertEqual(config["satwater"], 0.6)
        self.assertEqual(config["soil_moisture"], 0.5)


if __name__ == "__main__":
    unittest.main()
