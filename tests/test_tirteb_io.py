import csv
import math
import warnings
from pathlib import Path

import pytest

import tirteb_run
from base import canopy_closure, tirteb_radiation, tirteb_tirt
from base.io import read_defined, read_input

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "input_tirteb.csv"


def _config(**overrides):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        config = read_input(TEMPLATE)
    config["_input_base"] = str(ROOT)
    config.update({key: str(value) for key, value in overrides.items()})
    return config


def _portable_input(path, extra_rows=()):
    """Copy of the template with absolute file references, so it can live anywhere."""
    rows = []
    with TEMPLATE.open(newline="", encoding="utf-8-sig") as file:
        for row in csv.reader(file):
            if len(row) >= 3 and row[1].strip().lower() in tirteb_run.INPUT_FILE_KEYS:
                row[2] = str((ROOT / row[2]).resolve())
            rows.append(row)
    rows.extend(extra_rows)
    with path.open("w", newline="", encoding="utf-8") as file:
        csv.writer(file).writerows(rows)
    return path


def _run(input_path, output_path, **overrides):
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*TIRTEB parameter")
        result = tirteb_run.run(input_path, output_path, start_time=214.5, end_time=214.53,
                                observation_time_mode="all", output_intermediate="no", **overrides)
    with Path(result).open(newline="", encoding="utf-8") as file:
        return list(csv.DictReader(file))


# 1. Input snapshot never overwrites the input file.
def test_snapshot_beside_input_does_not_rewrite_input(tmp_path):
    source = _portable_input(tmp_path / "input.csv")
    original = source.read_bytes()
    _run(source, tmp_path / "output.csv")
    used = tmp_path / "input_used.csv"
    first = used.read_bytes()
    _run(source, tmp_path / "output.csv")
    assert source.read_bytes() == original
    assert used.read_bytes() == first
    assert used.read_text(encoding="utf-8").startswith(tirteb_run.SNAPSHOT_MARKER)
    assert first.count(b"[resolved_inputs]") == 1


def test_snapshot_keeps_foreign_input_csv_and_reuses_own_snapshot(tmp_path):
    source = _portable_input(tmp_path / "case.csv")
    results = tmp_path / "results"
    results.mkdir()
    foreign = results / "input.csv"
    foreign.write_text("user file\n", encoding="utf-8")
    _run(source, results / "output.csv")
    assert foreign.read_text(encoding="utf-8") == "user file\n"
    assert (results / "input_used.csv").is_file()

    other = tmp_path / "other"
    _run(source, other / "output.csv")
    snapshot = other / "input.csv"
    assert tirteb_run._is_input_snapshot(snapshot)
    size = snapshot.stat().st_size
    _run(source, other / "output.csv")
    assert snapshot.stat().st_size == size
    _run(snapshot, other / "output.csv")
    assert snapshot.stat().st_size == size
    assert (other / "input_used.csv").is_file()


# 2. Crown parameters follow TiRT: crown_rad_a, crown_rad_b, crown_height (base height).
def test_crown_keys_reach_every_consumer():
    settings = read_defined(None, _config(vegetation_model="crown", crown_rad_a=1.5, crown_rad_b=2.5,
                                          crown_height=1.2, crown_density=0.05))
    assert (settings["width"], settings["height"], settings["stand_density"]) == (1.5, 2.5, 0.05)
    assert settings["crown_base_height"] == 1.2
    assert tirteb_radiation.canopy(settings)[-1] == 1.2
    assert canopy_closure._crown(settings)[-1] == 1.2
    _, _, config = tirteb_tirt.tirt_config(settings, {"Tss": 300.0, "Tsh": 300.0, "Tls": 300.0, "Tlh": 300.0,
                                                      "sza": 30.0, "saa": 0.0})
    assert (config["crown_rad_a"], config["crown_rad_b"], config["crown_height"]) == (1.5, 2.5, 1.2)


def test_default_crown_geometry_is_unchanged():
    settings = read_defined(None, _config(vegetation_model="crown"))
    assert (settings["width"], settings["height"], settings["stand_density"], settings["crown_base_height"]) == (
        3.0, 3.0, 0.1, 3.0)


@pytest.mark.parametrize("key", ["crown_radius", "crown_width"])
def test_old_crown_keys_are_rejected(key):
    with pytest.raises(ValueError, match="crown_rad_a"):
        read_defined(None, _config(vegetation_model="crown", **{key: 2.0}))


# 3. TiRT slope names are aliases; unsupported TiRT geometry keys are rejected.
def test_slope_aliases():
    config = _config(surface_model="terrain", terrain_type=1)
    del config["terrain_slope"], config["terrain_aspect"]
    settings = read_defined(None, dict(config, slope_angle="20", slope_aspect="135"))
    assert (settings["terrain_slope"], settings["terrain_aspect"]) == (20.0, 135.0)
    same = read_defined(None, _config(slope_angle=30, terrain_slope=30))
    assert same["terrain_slope"] == 30.0
    with pytest.raises(ValueError, match="slope_angle"):
        read_defined(None, _config(slope_angle=20))


@pytest.mark.parametrize("key", ["geometry_mode", "geometry_file", "principal_vza", "principal_vza_step",
                                 "geometry_pair_mode"])
def test_unsupported_geometry_keys(key):
    with pytest.raises(ValueError, match=key):
        read_defined(None, _config(**{key: "0"}))


# 4. Unknown or unused keys warn with their name; the template is warning-free.
def test_template_is_warning_free():
    _config()


@pytest.mark.parametrize("key", ["lia", "aero_zo", "lidfa", "g", "doy", "soil_temperature_c",
                                 "urban_temperature_c", "urban_street_radiation_factor",
                                 "soil_depth_step", "urban_depth_step", "width", "stand_density"])
def test_unknown_keys_warn(tmp_path, key):
    path = _portable_input(tmp_path / "input.csv", [["extra", key, "1", ""]])
    with pytest.warns(UserWarning, match=f"'{key}'"):
        read_input(path)


# 6. Observation hours are clock hours.
def test_observation_hours():
    fractions = tirteb_run._observation_hour_fractions({"observation_hours": "1;12.5;24"})
    assert fractions == pytest.approx((1 / 24, 12.5 / 24, 0.0))
    assert tirteb_run._matches_observation_hour(214.0 + 1 / 24, fractions)
    assert not tirteb_run._matches_observation_hour(214.0, (12.5 / 24,))
    rounded_meteo_time = 214.5139  # 12:20 written with four decimals
    assert tirteb_run._matches_observation_hour(rounded_meteo_time, tirteb_run._observation_hour_fractions(
        {"observation_hours": "12.3333"}))
    assert not tirteb_run._matches_observation_hour(rounded_meteo_time, (12.5 / 24,))


@pytest.mark.parametrize("value", ["*.5", "25", "-1", "noon"])
def test_invalid_observation_hours(value):
    with pytest.raises(ValueError, match="observation_hours"):
        tirteb_run._observation_hour_fractions({"observation_hours": value})


# 7. Output columns: green band name, sunlit shares; absent leaves stay empty.
def test_output_columns(tmp_path):
    rows = _run(TEMPLATE, tmp_path / "output.csv")
    assert "green" in rows[0] and "red" not in rows[0]
    assert list(rows[0]) == list(tirteb_run.OUTPUT_FIELDS)
    for key in tirteb_run.BOOLEAN_DIAGNOSTICS:
        assert {row[key] for row in rows} <= {"0", "1"}
    for row in rows:
        assert 0.0 <= float(row["frac_soil_sunlit"]) <= 1.0
        assert 0.0 <= float(row["frac_leaf_sunlit"]) <= 1.0
        assert row["frac_roof_sunlit"] == "" or math.isnan(float(row["frac_roof_sunlit"]))
        assert math.isfinite(float(row["Tls_C"])) and math.isfinite(float(row["Tlh_C"]))


def test_urban_bare_keeps_both_components_but_no_leaves(tmp_path):
    rows = _run(TEMPLATE, tmp_path / "output.csv", surface_model="urban", vegetation_model="bare")
    for row in rows:
        for key in ("Troof_sunlit_C", "Troof_shaded_C", "Twall_sunlit_C", "Twall_shaded_C",
                    "Tstreet_sunlit_C", "Tstreet_shaded_C"):
            assert math.isfinite(float(row[key]))
        for key in ("frac_roof_sunlit", "frac_wall_sunlit", "frac_street_sunlit"):
            assert 0.0 <= float(row[key]) <= 1.0
        for key in ("Tls_C", "Tlh_C", "Tleaf_sunlit_C", "Tleaf_shaded_C", "frac_leaf_sunlit"):
            assert math.isnan(float(row[key]))
