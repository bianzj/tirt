"""Regression checks for wavelength units, conduction and EB result provenance."""

import csv
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from base.brightness import _urban_spectral_value, _vegetation_spectral_value
from base.eb import CP, MH20, MAIR, RAD_THRESHOLD, RHOA, SIGMASB, latent_heat, saturation_vapor_pressure
from base.evapo import _depth_differential
from base.io import read_defined, read_input, read_meteo
from base.pipeline import run_node
from base.structs import InputStructure
from tirteb_run import _write_input_snapshot, run


class SpectralAndConductionTests(unittest.TestCase):
    def test_urban_dynamic_observation_uses_resolved_sky_partition(self):
        from tirteb_run import _band_observation
        settings = {"surface_model": "urban", "shortwave_diffuse_fraction": -1}
        state = {"sza": 30, "saa": 0, "shortwave_diffuse_fraction": .32}
        with patch('tirteb_run._spectral_values', return_value=(.1, .2, .3)), \
             patch('tirteb_run.urban_directional_reflectance', return_value=.2) as urban, \
             patch('tirteb_run.component_proportions', return_value=np.array([1, 0, 0, 0])), \
             patch('tirteb_run.urban_component_fractions', return_value=np.ones(6)/6), \
             patch('tirteb_run.urban_canopy_fraction', return_value=.3):
            _band_observation(settings, {'shortwave': 700}, state, 'urban', .86, 20, 90)
        self.assertEqual(urban.call_args.args[0]['shortwave_diffuse_fraction'], .32)
        self.assertEqual(settings['shortwave_diffuse_fraction'], -1)

    def test_urban_query_converts_um_to_nm_for_optical_and_thermal_bands(self):
        settings = {
            "wavelengths": [400.0, 550.0, 860.0, 10500.0],
            "urban_reflectance_roof_spectrum": [0.08, 0.11, 0.16, 0.05],
            "urban_emissivity_roof_spectrum": [0.91, 0.92, 0.93, 0.95],
        }
        for wavelength, expected in ((0.55, 0.11), (0.86, 0.16)):
            self.assertAlmostEqual(_urban_spectral_value(
                settings, "urban_reflectance_roof_spectrum", wavelength, 0.0
            ), expected)
        self.assertAlmostEqual(_urban_spectral_value(
            settings, "urban_emissivity_roof_spectrum", 10.5, 0.0
        ), 0.95)

    def test_vegetation_optical_spectrum_uses_optical_part_of_full_grid(self):
        spectrum = np.linspace(0.0, 1.0, 2001)
        settings = {
            "wavelengths": list(range(400, 2401)) + [2500, 10500],
            "leaf_reflectance_spectrum": spectrum,
            "leaf_reflectance": 0.015,
        }
        self.assertAlmostEqual(_vegetation_spectral_value(
            settings, "leaf_reflectance_spectrum", "leaf_reflectance", 0.86
        ), spectrum[460])

    def test_depth_conduction_uses_subsurface_gradient_at_steady_surface(self):
        settings = {"soil_lambda": 1.55, "soil_depth_step": 0.05}
        # A steady surface still conducts heat into a colder first soil layer.
        state = {"depth_profile_sunlit": [300.0, 290.0, 280.0, 280.0]}
        self.assertAlmostEqual(_depth_differential(
            settings, state, 300.0, "depth_profile_sunlit"
        ), 310.0)
        state["depth_profile_sunlit"][0] = 285.0
        self.assertAlmostEqual(_depth_differential(
            settings, state, 300.0, "depth_profile_sunlit"
        ), 310.0)

    def test_snapshot_preserves_source_overrides_and_freezes_local_references(self):
        with tempfile.TemporaryDirectory(prefix="tirteb_snapshot_") as temp:
            root = Path(temp)
            source = root / "source"
            source.mkdir()
            meteo = source / "meteo.txt"
            meteo.write_text("1\n214.0139 20 15 900 1 0 350\n", encoding="utf-8")
            original_text = "# custom study\n[case],,,\ncase,meteo,meteo.txt,driver\nrun,lai,1.5,original\n"
            input_path = source / "input.csv"
            input_path.write_text(original_text, encoding="utf-8")
            config = read_input(input_path)
            config.update(lai="3.0", direction_source="direct", spectral_source="direct",
                          case_result_dir=str(root / "previous_case"))
            input_data = InputStructure.from_config(config, source)
            snapshot = root / "case" / "input.csv"
            _write_input_snapshot(input_path, snapshot, config,
                                  {"lai": 3.0, "case_result_dir": config["case_result_dir"]}, input_data)
            frozen_config = read_input(snapshot)
            self.assertEqual(snapshot.with_name("source_input.csv").read_text(encoding="utf-8"), original_text)
            self.assertEqual(frozen_config["lai"], "3.0")
            self.assertEqual(InputStructure.from_config(frozen_config, snapshot.parent).case_dir, snapshot.parent)
            self.assertFalse(Path(frozen_config["meteo"]).is_absolute())
            frozen_meteo = snapshot.parent / frozen_config["meteo"]
            expected = frozen_meteo.read_bytes()
            meteo.write_text("changed source", encoding="utf-8")
            self.assertEqual(frozen_meteo.read_bytes(), expected)
            # Re-running the snapshot itself neither copies a file onto itself
            # nor replaces the original source with the expanded snapshot.
            count = len(list(snapshot.parent.joinpath("inputs").iterdir()))
            _write_input_snapshot(snapshot, snapshot, frozen_config, {},
                                  InputStructure.from_config(frozen_config, snapshot.parent))
            self.assertEqual(len(list(snapshot.parent.joinpath("inputs").iterdir())), count)
            self.assertEqual(snapshot.with_name("source_input.csv").read_text(encoding="utf-8"), original_text)


class EnergyBalanceIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = read_input(ROOT / "input_tirteb.csv")
        cls.config.update(_input_base=str(ROOT), surface_model="plane",
                          temperature_method="robust", soil_temperature_method="direct")
        cls.settings = read_defined(None, cls.config)
        cls.records = read_meteo(ROOT / "data/defined/meteo.txt")
        cls.meteo = cls.records[0]

    def node(self, settings, max_iter):
        return run_node(settings, self.meteo, 0, 2024, 172, 39.9, 116.4,
                        max_iter=max_iter, utc_offset_hours=8)

    def test_fluxes_and_net_radiation_match_returned_temperature(self):
        settings = dict(self.settings)
        settings["soil_temperature_method"] = "robust"
        for limit in (1, 50):
            with self.subTest(max_iter=limit):
                state = self.node(settings, limit)
                temperature = state["Tss"]
                air_temperature = self.meteo["air_temperature_c"] + 273.15
                resistance = state["raa_soil"]
                expected_h = RHOA * CP / resistance * (temperature - air_temperature)
                expected_le = (
                    RHOA * latent_heat(temperature)
                    / (resistance + state["rss"])
                    * (MH20 / MAIR / self.meteo["pressure"])
                    * (saturation_vapor_pressure(temperature - 273.15) - self.meteo["vapor_pressure"])
                )
                expected_n = state["radiation_input_soil_sunlit"] - SIGMASB * temperature ** 4
                self.assertAlmostEqual(state["H_soilsunlit"], expected_h, places=10)
                self.assertAlmostEqual(state["LE_soilsunlit"], expected_le, places=10)
                self.assertAlmostEqual(state["N_soilsunlit"], expected_n, places=10)
                self.assertTrue(math.isfinite(state["energy_balance_residual_W_m2"]))
                if state["closed"]:
                    self.assertLessEqual(state["energy_balance_residual_W_m2"], RAD_THRESHOLD)
                if limit == 1:
                    self.assertFalse(state["closed"])
                    self.assertEqual(state["iterations"], 1)
                    self.assertGreater(state["energy_balance_residual_W_m2"], RAD_THRESHOLD)

    def test_soil_algorithm_override_does_not_replace_leaf_algorithm(self):
        self.assertEqual(self.settings["temperature_method"], "robust")
        self.assertEqual(self.settings["soil_temperature_method"], "direct")
        self.assertEqual(InputStructure.from_config(self.config, ROOT).soil_temperature_method, "direct")
        fallback = {**self.config, "temperature_method": "direct", "soil_temperature_method": ""}
        self.assertEqual(InputStructure.from_config(fallback, ROOT).soil_temperature_method, "direct")
        direct = self.node(dict(self.settings), 1)
        robust_settings = dict(self.settings)
        robust_settings["soil_temperature_method"] = "robust"
        robust = self.node(robust_settings, 1)
        self.assertNotAlmostEqual(direct["Tss"], robust["Tss"], places=6)
        self.assertAlmostEqual(direct["Tls"], robust["Tls"], places=10)
        self.assertAlmostEqual(direct["Tlh"], robust["Tlh"], places=10)

    def test_unclosed_status_survives_without_intermediate_csv(self):
        with tempfile.TemporaryDirectory(prefix="tirteb_physics_") as temp:
            output = Path(temp) / "output.csv"
            with patch("tirteb_run.read_defined", return_value=self.settings), self.assertWarnsRegex(RuntimeWarning, "能量平衡未闭合"):
                run(ROOT / "input_tirteb.csv", output,
                    start_time=self.meteo["time"], end_time=self.meteo["time"],
                    observation_time_mode="all", output_intermediate="no", max_iter=1)
            with output.open(newline="", encoding="utf-8") as file:
                row = next(csv.DictReader(file))
            self.assertEqual(row["closed"], "False")
            self.assertEqual(row["iterations"], "1")
            self.assertGreater(float(row["energy_balance_residual_W_m2"]), RAD_THRESHOLD)
            self.assertFalse(output.with_name("intermediate.csv").exists())

    def test_unmatched_observation_time_fails_before_output_changes(self):
        with tempfile.TemporaryDirectory(prefix="tirteb_physics_") as temp:
            output = Path(temp) / "output.csv"
            output.write_text("previous result\n", encoding="utf-8")
            with patch("tirteb_run.read_defined", return_value=self.settings), self.assertRaisesRegex(ValueError, "可选实际时间") as raised:
                run(ROOT / "input_tirteb.csv", output,
                    start_time=214.0, end_time=215.0,
                    observation_time_mode="time", observation_times="172.5")
            self.assertIn(str(self.meteo["time"]), str(raised.exception))
            self.assertEqual(output.read_text(encoding="utf-8"), "previous result\n")

    def test_default_observation_time_exists_in_default_simulation_window(self):
        targets = {float(item) for item in self.config["observation_times"].split(";")}
        start = float(self.config["start_time"])
        end = float(self.config["end_time"])
        self.assertTrue(any(start <= row["time"] <= end and row["time"] in targets
                            for row in self.records))


if __name__ == "__main__":
    unittest.main()
