import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from base.evapo import SIGMASB, _urban_component_incoming, _urban_longwave_inputs


class UrbanWallRadiationTests(unittest.TestCase):
    def test_wall_longwave_view_factors_are_closed(self):
        temperature = 300.0
        sky = SIGMASB * temperature ** 4
        settings = {
            "urban_building_length": 10.0,
            "urban_building_width": 10.0,
            "urban_building_height": 20.0,
            "urban_emissivity_roof": 1.0,
            "urban_emissivity_wall": 1.0,
            "urban_emissivity_street": 1.0,
            "leaf_emissivity": 1.0,
        }
        state = {
            "Troof_sunlit": temperature,
            "Troof_shaded": temperature,
            "Twall_sunlit": temperature,
            "Twall_shaded": temperature,
            "Tss": temperature,
            "Tsh": temperature,
        }

        result = _urban_longwave_inputs(settings, state, sky)

        self.assertAlmostEqual(result["LW_down_wall_sunlit"], sky, places=9)
        self.assertAlmostEqual(result["wall_sunlit"], sky, places=9)

    def test_sunlit_wall_projection_uses_illuminated_facades(self):
        settings = {
            "urban_building_length": 10.0,
            "urban_building_width": 10.0,
            "urban_building_height": 20.0,
            "urban_absorption_roof": 1.0,
            "urban_absorption_wall": 1.0,
            "urban_absorption_street": 1.0,
        }
        state = {
            "sza": 45.0,
            "saa": 45.0,
            "Troof_sunlit": 300.0,
            "Troof_shaded": 300.0,
            "Twall_sunlit": 300.0,
            "Twall_shaded": 300.0,
            "Tss": 300.0,
            "Tsh": 300.0,
        }
        netrad = {
            "directVrad_incident": 500.0,
            "diffuseVrad_soil": 0.0,
            "longwave_down": 0.0,
        }

        _urban_component_incoming(settings, state, netrad)

        expected_projection = math.sqrt(0.5)
        self.assertAlmostEqual(state["wall_orientation"], expected_projection, places=9)
        self.assertAlmostEqual(state["wall_direct_factor"], expected_projection, places=9)
        self.assertAlmostEqual(state["wall_sunlit_fraction"], 0.5, places=9)
        self.assertAlmostEqual(
            state["SW_direct_wall_sunlit"], 500.0 * expected_projection, places=9
        )


if __name__ == "__main__":
    unittest.main()
