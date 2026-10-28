"""Soil water uses one scale: volume fraction (m3 m-3) with 0.25 as BSM reference."""
from pathlib import Path

import pytest

from base.io import _read_optical_coefficients, read_defined, read_input
from base.soilopt import bsm, bsm_water_mu, soil_tir_emissivity, volumetric_water

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "input_tirteb.csv"


def _config(**overrides):
    config = read_input(TEMPLATE)
    config["_input_base"] = str(ROOT)
    config.update({key: str(value) for key, value in overrides.items()})
    return config


@pytest.mark.parametrize("value, expected", [
    (0.25, 0.25), (25, 0.25), (0.45, 0.45), (45, 0.45), (1.0, 1.0),
])
def test_volumetric_water(value, expected):
    assert volumetric_water(value) == pytest.approx(expected)


def test_bsm_water_mu_is_referenced_to_0_25():
    assert bsm_water_mu(0.05) == pytest.approx(0.0)
    assert bsm_water_mu(0.25) == pytest.approx(0.8)
    assert bsm_water_mu(0.45) == pytest.approx(1.6)
    # Legacy percent input must land on the same mu.
    assert bsm_water_mu(25) == pytest.approx(bsm_water_mu(0.25))


def test_bsm_accepts_fraction_and_percent_equally():
    coeff = _read_optical_coefficients(ROOT / "data" / "defined" / "optipar_new.csv")
    soil_optical = [0.25, 25.0, 45.0]
    dry = bsm(coeff, soil_optical, 0.05)
    at_reference = bsm(coeff, soil_optical, 0.25)
    at_percent = bsm(coeff, soil_optical, 25.0)
    wet = bsm(coeff, soil_optical, 0.45)
    for index in (150, 460, 1800):  # 550, 860 and 2200 nm
        assert at_reference[index] == pytest.approx(at_percent[index])
        assert dry[index] > at_reference[index] > wet[index]


def test_tir_water_correction_accepts_percent_references():
    assert soil_tir_emissivity(
        25, enabled=1, dry_moisture=5, wet_moisture=45,
    ) == pytest.approx(soil_tir_emissivity(
        0.25, enabled=1, dry_moisture=0.05, wet_moisture=0.45,
    ))


def test_defined_structure_normalizes_every_soil_water_input():
    settings = read_defined(None, _config(
        soil_moisture=25, satwater=45,
        soil_tir_moisture_dry=5, soil_tir_moisture_wet=45,
    ))
    assert settings["soil_moisture"] == pytest.approx(0.25)
    assert settings["satwater"] == pytest.approx(0.45)
    assert settings["soil_tir_moisture_dry"] == pytest.approx(0.05)
    assert settings["soil_tir_moisture_wet"] == pytest.approx(0.45)
