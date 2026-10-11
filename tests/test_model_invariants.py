"""Physical invariants of the directional models (isothermal blackbody, finiteness, emissivity response)."""

import itertools
import sys
import unittest
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import run as R  # noqa: E402
from base.inputio import read_input  # noqa: E402

SURFACES = ("plane", "slope", "terrain", "urban")
VEGETATIONS = ("bare", "hom", "row", "crown")
WAVELENGTHS = np.array([8.0, 10.5, 12.0])


def _geometries(config, surface):
    for vza, sza, raa in itertools.product(range(0, 91, 15), range(0, 86, 15), (0, 90, 180)):
        geometry = (float(vza), float(raa), float(sza), 0.0)
        if surface == "slope" and R._slope_backface(config, geometry):
            continue
        yield geometry


class ModelInvariantTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore", RuntimeWarning)
        cls.config = read_input(R.DEFAULT_INPUT)
        cls.spectra = R._spectra(cls.config, R.DEFAULT_INPUT, WAVELENGTHS)

    def _results(self, config, spectra, surface, vegetation):
        for geometry in _geometries(config, surface):
            for band, wavelength in enumerate(WAVELENGTHS):
                band_spectra = {key: values[band] for key, values in spectra.items()}
                yield R._surface_result(config, surface, vegetation, float(wavelength), geometry, band_spectra)

    def test_all_models_are_finite(self):
        for surface, vegetation in itertools.product(SURFACES, VEGETATIONS):
            with self.subTest(model=f"{surface}+{vegetation}"):
                for radiance, temperature in self._results(self.config, self.spectra, surface, vegetation):
                    self.assertTrue(np.isfinite(radiance) and np.isfinite(temperature))

    def test_isothermal_blackbody_returns_its_temperature(self):
        config = dict(self.config, hapke_enable="0")
        for key in config:
            if key.endswith("_c") and ("sunlit" in key or "shaded" in key):
                config[key] = "26.85"
        black = {key: np.ones_like(values) for key, values in self.spectra.items()}
        for surface, vegetation in itertools.product(SURFACES, VEGETATIONS):
            with self.subTest(model=f"{surface}+{vegetation}"):
                for _, temperature in self._results(config, black, surface, vegetation):
                    self.assertAlmostEqual(temperature, 300.0, places=6)

    def test_terrain_vegetation_responds_to_emissivity(self):
        geometry = (30.0, 0.0, 30.0, 0.0)
        band = {key: values[1] for key, values in self.spectra.items()}
        low = dict(band, leaf_e=0.9, terrain_e=0.9)
        warm, _ = R._surface_result(self.config, "terrain", "hom", 10.5, geometry, band)
        cool, _ = R._surface_result(self.config, "terrain", "hom", 10.5, geometry, low)
        self.assertLess(cool, warm)

    def test_flat_bare_soil_has_no_shadows(self):
        band = {key: values[1] for key, values in self.spectra.items()}
        config = dict(self.config, hapke_enable="0")
        soil_sun, soil_shade, _, _ = R._vegetation_temperatures(config)
        sunlit = band["soil_e"] * R.planck(10.5, soil_sun)
        shaded = band["soil_e"] * R.planck(10.5, soil_shade)
        for vza, vaa in ((0.0, 0.0), (45.0, 180.0), (75.0, 90.0)):
            radiance, _ = R._surface_result(config, "plane", "bare", 10.5, (vza, vaa, 30.0, 0.0), band)
            self.assertAlmostEqual(radiance, sunlit, places=9)
        radiance, _ = R._surface_result(dict(config, lai="0"), "plane", "hom", 10.5, (30.0, 0.0, 30.0, 0.0), band)
        self.assertAlmostEqual(radiance, sunlit, places=9)
        radiance, _ = R._surface_result(config, "plane", "bare", 10.5, (30.0, 0.0, 95.0, 0.0), band)
        self.assertAlmostEqual(radiance, shaded, places=9)
        # Sun behind a north-facing 30 degree slope (sun from the south, 70 degree zenith).
        slope = dict(config, slope_angle="30", slope_aspect="0")
        radiance, _ = R._surface_result(slope, "slope", "bare", 10.5, (10.0, 0.0, 70.0, 180.0), band)
        self.assertAlmostEqual(radiance, shaded, places=9)

    def test_row_fractions_match_cross_section_integral(self):
        from base.proportion import proportion_bidirectional_row_one

        k, width, gap, height, sza = 1.5, 0.5, 0.5, 1.0, 30.0
        period = width + gap
        x0 = (np.arange(2000) + 0.5) / 2000 * period
        z = np.linspace(0.0, height, 2001)

        def tau(slope, cos_zenith):
            position = x0[:, None] + slope * z[None, :]
            return k / cos_zenith * np.mean((np.mod(position, period) < width), axis=1) * height

        def row(vza, vaa, row_az):
            one = lambda value: np.array([float(value)])
            parts = proportion_bidirectional_row_one(1.5, 1.0e-6, width, gap, height, one(vza), one(vaa), one(sza), one(0.0), one(row_az))
            return [float(np.asarray(part).reshape(-1)[0]) for part in parts]

        # A vanishing hspot (leaf size) makes the leaf-scale overlap weight 0, leaving the exact cross-section overlap.
        rd = np.deg2rad
        for row_az in (0.0, 90.0):
            for vza, vaa in ((0.0, 0.0), (45.0, 180.0), (60.0, 180.0), (60.0, 90.0), (75.0, 270.0)):
                tv = tau(np.tan(rd(vza)) * np.sin(rd(vaa - row_az)), np.cos(rd(vza)))
                ts = tau(np.tan(rd(sza)) * np.sin(rd(0.0 - row_az)), np.cos(rd(sza)))
                soil_sunlit, soil_shaded, leaf_sunlit, leaf_shaded = row(vza, vaa, row_az)
                with self.subTest(row_az=row_az, vza=vza, vaa=vaa):
                    self.assertAlmostEqual(soil_sunlit, float(np.mean(np.exp(-tv - ts))), delta=0.01)
                    self.assertAlmostEqual(soil_sunlit + soil_shaded, float(np.mean(np.exp(-tv))), delta=0.01)
                    self.assertAlmostEqual(soil_sunlit + soil_shaded + leaf_sunlit + leaf_shaded, 1.0, places=9)
                    self.assertGreaterEqual(min(soil_sunlit, soil_shaded, leaf_sunlit, leaf_shaded), -1.0e-12)

    def test_crown_fractions_match_poisson_crown_integral(self):
        from base.proportion import proportion_bidirectional_crown_one

        lai, std, sza = 1.5, 0.02, 30.0
        rng = np.random.default_rng(0)

        def chord(offsets, zenith, azimuth, a, b):
            # Path length through a spheroid centred at `offsets` (crown centre minus ground point) of the ray leaving the ground.
            d = np.array([np.sin(zenith) * np.cos(azimuth), np.sin(zenith) * np.sin(azimuth), np.cos(zenith)])
            scale = np.array([1 / a, 1 / a, 1 / b])
            p, q = -offsets * scale, d * scale
            A, B, C = q @ q, p @ q, np.sum(p * p, axis=1) - 1
            return 2 * np.sqrt(np.maximum(B * B - A * C, 0.0)) / A

        def crown(a, b, base, vza, vaa):
            one = lambda value: np.array([float(value)])
            parts = proportion_bidirectional_crown_one(lai, std, 1.0e-6, b, a, one(vza), one(sza), one(vaa), hcb=base)
            return [float(np.asarray(part).reshape(-1)[0]) for part in parts]

        # A vanishing hspot removes the leaf-scale overlap; soil gaps are then exp(-std * integral of (1 - T)) over crown positions.
        rd = np.deg2rad
        for a, b, base in ((3.0, 3.0, 3.0), (1.5, 4.0, 1.0)):
            extinction = 0.5 * lai / (std * 4 / 3 * np.pi * a * a * b)
            for vza, vaa in ((0.0, 0.0), (20.0, 0.0), (45.0, 180.0), (60.0, 90.0)):
                half = (base + 2 * b) * np.tan(rd(max(vza, sza))) + a + 0.5
                offsets = np.column_stack([rng.uniform(-half, half, 400000), rng.uniform(-half, half, 400000), np.full(400000, base + b)])
                tau_v = extinction * chord(offsets, rd(vza), rd(vaa), a, b)
                tau_s = extinction * chord(offsets, rd(sza), 0.0, a, b)
                area = (2 * half) ** 2
                soil = np.exp(-std * area * np.mean(-np.expm1(-tau_v)))
                soil_sunlit = np.exp(-std * area * np.mean(-np.expm1(-tau_v - tau_s)))
                fractions = crown(a, b, base, vza, vaa)
                with self.subTest(crown=(a, b, base), vza=vza, vaa=vaa):
                    self.assertAlmostEqual(fractions[0] + fractions[1], soil, delta=0.006)
                    self.assertAlmostEqual(fractions[0], soil_sunlit, delta=0.006)
                    self.assertAlmostEqual(sum(fractions), 1.0, places=9)
                    self.assertGreaterEqual(min(fractions), -1.0e-12)
        # At the hotspot nothing visible is shaded.
        _, soil_shaded, _, leaf_shaded = crown(3.0, 3.0, 3.0, sza, 0.0)
        self.assertAlmostEqual(soil_shaded, 0.0, places=9)
        self.assertAlmostEqual(leaf_shaded, 0.0, places=9)


if __name__ == "__main__":
    unittest.main()
