"""TIRTEB thermal output is TiRT: the same input file and component temperatures give the same radiance."""

import csv
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import run as tirt  # noqa: E402
import tirteb_run as tirteb  # noqa: E402
from base.inputio import read_input  # noqa: E402
from base.io import read_defined  # noqa: E402

INPUT = ROOT / "input_tirteb.csv"
TEMPERATURES_C = {
    "soil": (41.0, 27.0), "leaf": (31.0, 28.0), "roof": (46.0, 29.0), "wall": (39.0, 30.0), "street": (44.0, 28.5),
}
STATE_KEYS = {
    "soil": ("Tss", "Tsh"), "leaf": ("Tls", "Tlh"), "roof": ("Troof_sunlit", "Troof_shaded"),
    "wall": ("Twall_sunlit", "Twall_shaded"), "street": ("Tstreet_sunlit", "Tstreet_shaded"),
}
SUN = (35.0, 140.0)
CASES = (
    ("plane", "bare"), ("plane", "hom"), ("plane", "row"), ("plane", "crown"),
    ("slope", "hom"), ("slope", "crown"), ("terrain", "bare"), ("terrain", "crown"),
    ("urban", "bare"), ("urban", "crown"),
)


class TirtebUsesTirtTests(unittest.TestCase):
    def test_same_input_gives_same_radiance(self):
        warnings.simplefilter("ignore", RuntimeWarning)
        state = {"sza_global": SUN[0], "saa_global": SUN[1], "sunlit_surface": True}
        temperatures = {}
        for component, (sunlit, shaded) in TEMPERATURES_C.items():
            state[STATE_KEYS[component][0]] = sunlit + 273.15
            state[STATE_KEYS[component][1]] = shaded + 273.15
            temperatures[f"{component}_sunlit_c"] = sunlit
            temperatures[f"{component}_shaded_c"] = shaded
        temperatures["terrain_sunlit_c"], temperatures["terrain_shaded_c"] = TEMPERATURES_C["soil"]
        work = Path(tempfile.mkdtemp(prefix="tirteb-tirt-"))
        for surface, vegetation in CASES:
            overrides = {
                "surface_model": surface, "vegetation_model": vegetation,
                "spectral_source": "direct", "wavelengths": "8;10.5;12",
                "emissivity_leaf": "0.97;0.985;0.98", "emissivity_soil": "0.93;0.955;0.96",
                "geometry_mode": "0", "geometry_pair_mode": "zip",
                "vza": "0;30;45;60", "vaa": "0;140;320;230", "sza": SUN[0], "saa": SUN[1],
            }
            with self.subTest(model=f"{surface}+{vegetation}"):
                output = tirt.run(INPUT, work / f"{surface}_{vegetation}.csv", **overrides, **temperatures)
                with output.open(encoding="utf-8") as file:
                    rows = list(csv.DictReader(file))
                self.assertGreaterEqual(len(rows), 9)
                config = read_input(INPUT)
                config.update({key: str(item) for key, item in overrides.items()})
                settings = read_defined(config)
                for row in rows:
                    radiance, _ = tirteb._thermal_radiance(
                        config, settings, state, float(row["wavelength_um"]), float(row["vza"]), float(row["vaa"]))
                    self.assertAlmostEqual(radiance, float(row["radiance"]), places=9)


if __name__ == "__main__":
    unittest.main()
