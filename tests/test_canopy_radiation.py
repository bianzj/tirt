"""Energy conservation of the TIRTEB canopy radiation and of its scene fluxes."""

import csv
import math
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from base.radiation import flux_weights  # noqa: E402
from base.rt_transfer import _canopy_transfer, canopy_longwave, canopy_shortwave  # noqa: E402
from tirteb_run import run  # noqa: E402

SIGMA = 5.67e-8


class CanopyRadiationTests(unittest.TestCase):
    def test_isothermal_canopy_absorbs_its_own_emission(self):
        temperature = 300.0
        sky = SIGMA * temperature ** 4
        for leaf_e, soil_e in ((1.0, 1.0), (0.98, 0.94)):
            for lai, extinction in ((0.5, 0.0), (3.0, 0.6), (6.0, 1.2)):
                sunlit, shaded, soil = canopy_longwave(lai, leaf_e, soil_e, sky, temperature, temperature,
                                                       temperature, temperature, extinction, 0.3)
                with self.subTest(leaf_e=leaf_e, lai=lai):
                    self.assertAlmostEqual(sunlit / (2.0 * leaf_e * sky), 1.0, places=3)
                    self.assertAlmostEqual(shaded / (2.0 * leaf_e * sky), 1.0, places=3)
                    self.assertAlmostEqual(soil / (soil_e * sky), 1.0, places=3)

    def test_black_canopy_absorbs_what_it_does_not_transmit(self):
        esun, esky, lai, sza = 500.0, 200.0, 3.0, 30.0
        t = _canopy_transfer(lai, 0.0, 0.0, sza)
        leaf_diffuse, soil_absorbed, soil_diffuse = canopy_shortwave(lai, 0.0, 0.0, 0.0, sza, esun, esky)
        self.assertAlmostEqual(leaf_diffuse * lai, esky * (1.0 - t["tau_dd"]), places=6)
        self.assertAlmostEqual(soil_diffuse, esky * t["tau_dd"], places=6)
        self.assertAlmostEqual(soil_absorbed, soil_diffuse, places=6)

    def test_flux_weights_cover_ground_and_leaf_area(self):
        settings = {"lai": 3.0, "vegetation_model": "hom", "hspot": 0.05}
        for sza in (0.0, 30.0, 60.0, 85.0, 95.0):
            soil_sunlit, soil_shaded, lai_sunlit, lai_shaded = flux_weights(settings, sza, 0.0)
            with self.subTest(sza=sza):
                self.assertAlmostEqual(soil_sunlit + soil_shaded, 1.0, places=9)
                self.assertAlmostEqual(lai_sunlit + lai_shaded, 3.0, places=9)
                self.assertGreaterEqual(min(soil_sunlit, lai_sunlit), 0.0)
        self.assertAlmostEqual(flux_weights(settings, 30.0, 0.0)[2], (1.0 - math.exp(-0.5 * 3.0 / math.cos(math.radians(30.0)))) / (0.5 / math.cos(math.radians(30.0))), delta=0.05)

    def test_scene_net_radiation_is_partitioned(self):
        warnings.simplefilter("ignore", RuntimeWarning)
        output = run(ROOT / "input_tirteb.csv", Path(tempfile.mkdtemp(prefix="tirteb-eb-")) / "output.csv",
                     start_time=214.40, end_time=214.60, observation_time_mode="all", geometry_mode=0, vza=0, vaa=0,
                     spectral_source="direct", wavelengths="10.5")
        with output.open(encoding="utf-8") as file:
            rows = list(csv.DictReader(file))
        self.assertTrue(rows)
        for row in rows:
            residual = float(row["Rn"]) - float(row["LE"]) - float(row["H"]) - float(row["G"])
            with self.subTest(time=row["time"]):
                # Each component closes to 5 W m-2 per unit soil or leaf area; the scene has 1 + LAI of them.
                self.assertLess(abs(residual), 5.0 * 3.0)


if __name__ == "__main__":
    unittest.main()
