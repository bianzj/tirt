"""TIRTEB input/entry fixes of 2026-10-28: sky reflection, key table, compare summary, observation times."""
import csv
import warnings
from pathlib import Path

import pytest

import tirteb_run
from base import tirteb_tirt
from base.io import read_defined, read_input
from base.structs import InputStructure, _time_numbers

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "input_tirteb.csv"


def _config(**overrides):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        config = read_input(TEMPLATE)
    config["_input_base"] = str(ROOT)
    config.update({key: str(value) for key, value in overrides.items()})
    return config


# 1. sky_temperature_k / sky_temperature_c reach the TiRT observation.
def test_sky_temperature_is_passed_through():
    # Default: derived from the measured downwelling longwave at every step.
    assert read_defined(None, _config())["sky_temperature_k"] == "auto"
    assert read_defined(None, _config(sky_temperature_k=""))["sky_temperature_k"] == "auto"
    assert read_defined(None, _config(sky_temperature_k="none"))["sky_temperature_k"] is None
    assert read_defined(None, _config(sky_temperature_k=280))["sky_temperature_k"] == 280.0
    assert read_defined(None, _config(sky_temperature_c=10))["sky_temperature_k"] == pytest.approx(283.15)
    with pytest.raises(ValueError, match="sky_temperature"):
        read_defined(None, _config(sky_temperature_k=-5))


def test_auto_sky_uses_measured_longwave():
    state = {"Tss": 300.0, "Tsh": 300.0, "Tls": 300.0, "Tlh": 300.0, "sza": 30.0, "saa": 0.0}
    settings = tirteb_tirt.tir_band_settings(read_defined(None, _config()), 10.5)
    fixed = tirteb_tirt.tir_band_settings(read_defined(None, _config(sky_temperature_k=280.0)), 10.5)
    down = 5.67e-8 * 280.0 ** 4
    auto = tirteb_tirt.tirt_brightness(settings, dict(state, longwave_down=down), 10.5, 0.0, 0.0)
    assert auto == pytest.approx(tirteb_tirt.tirt_brightness(fixed, state, 10.5, 0.0, 0.0))


def test_sky_reflection_adds_radiance():
    state = {"Tss": 305.0, "Tsh": 298.0, "Tls": 300.0, "Tlh": 297.0, "sza": 30.0, "saa": 0.0}
    black = tirteb_tirt.tir_band_settings(read_defined(None, _config(sky_temperature_k="none")), 10.5)
    sky = tirteb_tirt.tir_band_settings(read_defined(None, _config(sky_temperature_k=290)), 10.5)
    black_radiance, black_bt = tirteb_tirt.tirt_brightness(black, state, 10.5, 0.0, 0.0)
    sky_radiance, sky_bt = tirteb_tirt.tirt_brightness(sky, state, 10.5, 0.0, 0.0)
    assert sky_radiance > black_radiance and sky_bt > black_bt


def test_soil_moisture_is_capped_at_saturation():
    settings = read_defined(None, _config(soil_moisture=0.55, satwater=0.45))
    assert settings["soil_moisture"] == pytest.approx(0.45)


def test_aerodynamic_height_is_an_input_key():
    settings = read_defined(None, _config(aerodynamic_height=7.5))
    assert settings["aerodynamic_height"] == pytest.approx(7.5)


# 2. Key table: read keys do not warn as unknown; header rows are not parameters.
@pytest.mark.parametrize("key", ["sky_temperature_k", "sky_temperature_c", "aero_cr", "urban_temperature_method"])
def test_consumed_keys_do_not_warn(tmp_path, key):
    path = tmp_path / "input.csv"
    path.write_text(f"section,parameter,value\nextra,{key},1\n", encoding="utf-8")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        values = read_input(path)
    assert "parameter" not in values and values[key] == "1"


@pytest.mark.parametrize("key", ["aero_zo", "lidfa", "doy", "urban_street_radiation_factor", "terrain_base_model"])
def test_parsed_but_unused_keys_say_so(tmp_path, key):
    path = tmp_path / "input.csv"
    path.write_text(f"extra,{key},1\n", encoding="utf-8")
    with pytest.warns(UserWarning, match=f"'{key}' is read but has no effect"):
        read_input(path)


# 3. Method comparison summary identifies each method combination.
def test_comparison_summary_has_method_columns(tmp_path, monkeypatch):
    def fake_run(input_path, **overrides):
        output = Path(overrides["case_result_dir"]) / "output.csv"
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=tirteb_run.OUTPUT_FIELDS)
            writer.writeheader()
            writer.writerow({"node": 0, "time": 214.5})
        return output

    monkeypatch.setattr(tirteb_run, "run", fake_run)
    summary = tirteb_run.run_comparison(TEMPLATE, tmp_path)
    with summary.open(newline="", encoding="utf-8") as file:
        rows = list(csv.DictReader(file))
    assert list(rows[0])[:3] == ["soil_heat_method", "stomatal_method", "soil_temperature_method"]
    combinations = {(row["soil_heat_method"], row["stomatal_method"], row["soil_temperature_method"]) for row in rows}
    assert len(rows) == len(combinations) == 12


# 4. Observation times: invalid items raise; matching uses the half-minute tolerance.
@pytest.mark.parametrize("value", ["214.5;noon", "*.5", "214.5,nan"])
def test_invalid_observation_times_raise(value):
    with pytest.raises(ValueError, match="observation_times"):
        _time_numbers({"observation_times": value}, "observation_times")
    with pytest.raises(ValueError, match="observation_times"):
        InputStructure.from_config({"observation_times": value}, ROOT)


def test_observation_time_matching_uses_tolerance():
    assert _time_numbers({"observation_times": "214.5; 215.25"}, "observation_times") == (214.5, 215.25)
    assert tirteb_run._matches_observation_time(214.51391, (214.5139,))
    assert tirteb_run._matches_observation_time(214.5139 + 0.4 / 1440, (214.5139,))
    assert not tirteb_run._matches_observation_time(214.5139 + 1.0 / 1440, (214.5139,))


# 5. Default crown density matches TiRT.
def test_default_crown_density_matches_tirt():
    config = _config(vegetation_model="crown")
    del config["crown_density"]
    settings = read_defined(None, config)
    assert settings["stand_density"] == settings["crown_density"] == 0.02
    assert InputStructure.from_config({}, ROOT).crown_density == 0.02
