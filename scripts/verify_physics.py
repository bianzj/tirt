"""Run physical and input regressions without requiring a browser runtime."""

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

PATTERNS = (
    "test_tirt_physics.py",
    "test_tirt_input_snapshot.py",
    "test_tirt_gui_inputs.py",
    "test_tirt_dynamic_gui.py",
    "test_tirteb_physics.py",
    "test_urban_wall_radiation.py",
    "test_optical.py",
    "test_optical_crown.py",
    "test_thermal_crown.py",
    "test_forest_terrain.py",
    "test_optical_urban_profile.py",
    "test_optical_rows_terrain_source.py",
    "test_optical_terrain_ports.py",
    "test_thermal_row.py",
    "test_optical_crown_intrinsic.py",
    "test_thermal_forest_exchange.py",
    "test_forest_unified.py",
    "test_urban_forest_hotspot.py",
)


def main():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite(
        loader.discover(str(ROOT / "tests"), pattern=pattern)
        for pattern in PATTERNS
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
