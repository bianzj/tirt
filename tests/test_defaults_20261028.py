"""Round-3 defaults: unified terrain meaning, aerodynamic canopy top, Cca, sky temperature."""
import warnings
from pathlib import Path

import pytest

import run
from base.io import read_defined, read_input

ROOT = Path(__file__).resolve().parents[1]


def _settings(**overrides):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        config = read_input(ROOT / "input_tirteb.csv")
    config["_input_base"] = str(ROOT)
    config.update({key: str(value) for key, value in overrides.items()})
    return read_defined(None, config)


def test_terrain_means_compound_terrain_in_both_models():
    # TIRTEB: terrain without terrain_type is compound terrain; slope is the single slope.
    config = _settings(surface_model="terrain")
    config_without_type = {k: v for k, v in config.items()}
    assert _settings(surface_model="terrain")["terrain_type"] == 2
    assert _settings(surface_model="slope")["terrain_type"] == 1
    # TiRT: the legacy alias terrain_type=1 now selects the slope for every band.
    assert run._canonical({"surface_model": "terrain"})[0] == "terrain"
    assert run._canonical({"surface_model": "terrain", "terrain_type": "1"})[0] == "slope"
    with pytest.raises(ValueError, match="terrain_type"):
        run._canonical({"surface_model": "terrain", "terrain_type": "3"})
    assert config_without_type


def test_aerodynamic_height_is_canopy_top():
    crown = _settings(vegetation_model="crown", crown_height=4, crown_rad_b=3)
    assert crown["aerodynamic_height"] == pytest.approx(4 + 2 * 3)
    assert crown["height"] == pytest.approx(3)  # radiative vertical semiaxis is unchanged
    row = _settings(vegetation_model="row", row_height=1.5)
    assert row["aerodynamic_height"] == pytest.approx(1.5)
    hom = _settings(vegetation_model="hom", height=2.5)
    assert hom["aerodynamic_height"] == pytest.approx(2.5)
    assert _settings(vegetation_model="crown", aerodynamic_height=12)["aerodynamic_height"] == 12.0


def test_carotenoids_default_lowers_visible_leaf_reflectance():
    default = _settings()
    without = _settings(cca=0)
    i = 100  # 500 nm on the 400 + i grid
    assert default["leaf_reflectance_spectrum"][i] < without["leaf_reflectance_spectrum"][i]
    j = 460  # 860 nm: carotenoids do not absorb in the NIR
    assert default["leaf_reflectance_spectrum"][j] == pytest.approx(without["leaf_reflectance_spectrum"][j], abs=1e-9)


def test_tirt_sky_default_and_black_option():
    assert run._sky_temperature({}) == run.DEFAULT_SKY_TEMPERATURE_K
    assert run._sky_temperature({"sky_temperature_k": "none"}) is None
    assert run._sky_temperature({"sky_temperature_k": "270"}) == 270.0
