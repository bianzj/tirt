import csv
import math
import shutil
import warnings
from pathlib import Path

import pytest

import run
from base import optical_urban, urban
from base.inputio import SNAPSHOT_MARKER, read_input
from base.optical import InvisibleDirection, reflectance

ROOT = Path(__file__).resolve().parents[1]
SINGLE_BAND = dict(wavelengths="0.67", leaf_reflectance="0.05", leaf_transmittance="0.03",
                   soil_reflectance="0.12", street_reflectance="0.12", wall_reflectance="0.2",
                   roof_reflectance="0.15")


def _rows(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _run_quiet(*args, **overrides):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        path = run.run(*args, **overrides)
    return path, [str(item.message) for item in caught]


def test_snapshot_does_not_overwrite_user_input_csv(tmp_path):
    case = tmp_path / "t1"
    case.mkdir()
    shutil.copy(ROOT / "input.csv", case / "my.csv")
    user_file = case / "input.csv"
    user_file.write_text("my own notes\n", encoding="utf-8")

    run.run(case / "my.csv", "out.csv", vza="0;30", vaa="0")

    assert user_file.read_text(encoding="utf-8") == "my own notes\n"
    snapshot = case / "input_used.csv"
    assert snapshot.read_text(encoding="utf-8").startswith(SNAPSHOT_MARKER)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        replay = read_input(snapshot)
    assert replay["surface_model"] == "plane" and replay["vza"] == "0;30"


def test_existing_snapshot_is_refreshed_in_place(tmp_path):
    shutil.copy(ROOT / "input.csv", tmp_path / "my.csv")
    run.run(tmp_path / "my.csv", tmp_path / "out" / "output.csv", vza="0", vaa="0")
    run.run(tmp_path / "my.csv", tmp_path / "out" / "output.csv", vza="30", vaa="0")
    assert read_input(tmp_path / "out" / "input.csv")["vza"] == "30"
    assert not (tmp_path / "out" / "input_used.csv").exists()


def test_invalid_optical_settings_are_not_reported_as_hidden_directions(tmp_path):
    with pytest.raises(ValueError, match="optical_terrain_method: typo"):
        run.run(ROOT / "input.csv", tmp_path / "output.csv", surface_model="terrain",
                vegetation_model="crown", canopy_closure="legacy", optical_terrain_method="TYPO",
                vza="0", vaa="0", **SINGLE_BAND)


def test_hidden_slope_direction_is_omitted(tmp_path):
    with pytest.warns(run.BackfaceOmittedWarning):
        path = run.run(ROOT / "input.csv", tmp_path / "output.csv", surface_model="slope",
                       vegetation_model="hom", slope_angle="30", slope_aspect="0",
                       vza="0;80", vaa="180", **SINGLE_BAND)
    rows = _rows(path)
    assert [float(row["vza"]) for row in rows] == [0.0]
    assert issubclass(InvisibleDirection, ValueError)
    settings = dict(lai=1.5, wavelengths=[670.0], leaf_reflectance=0.05, leaf_transmittance=0.03,
                    soil_reflectance=0.12, surface_model="plane", vegetation_model="hom")
    with pytest.raises(InvisibleDirection):
        reflectance(settings, 670.0, 90.0, 0.0, 30.0, 0.0)


def _urban_case(tmp_path, monkeypatch, density, height="20"):
    seen = {"optical": set(), "thermal": set()}
    building_shape = optical_urban.building_shape
    set_structural_input = urban.Urban.set_structural_input

    def optical_shape(settings):
        shape = building_shape(settings)
        seen["optical"].add((shape[2], shape[3]))
        return shape

    def thermal_shape(self, shapes):
        seen["thermal"].update((float(row[2]), float(row[3])) for row in shapes)
        return set_structural_input(self, shapes)

    monkeypatch.setattr(optical_urban, "building_shape", optical_shape)
    monkeypatch.setattr(urban.Urban, "set_structural_input", thermal_shape)
    path = run.run(ROOT / "input.csv", tmp_path / f"urban_{density}_{height}" / "output.csv",
                   surface_model="urban", vegetation_model="bare", building_density=density,
                   building_height=height, vza="30", vaa="0", wavelengths="0.67;10.5",
                   leaf_reflectance="0.05;0.015", leaf_transmittance="0.03;0.0",
                   soil_reflectance="0.12;0.045", street_reflectance="0.12;0.045",
                   wall_reflectance="0.2;0.07", roof_reflectance="0.15;0.05")
    brf = [float(row["reflectance_brf"]) for row in _rows(path) if float(row["wavelength_um"]) < 2.5]
    return brf[0], seen


def test_urban_building_parameters_reach_the_optical_model(tmp_path, monkeypatch):
    results = {density: _urban_case(tmp_path, monkeypatch, density) for density in ("0", "0.0015", "0.003")}
    values = [brf for brf, _ in results.values()]
    assert all(math.isfinite(item) for item in values)
    assert len({round(item, 12) for item in values}) == 3
    _, seen = results["0.003"]
    assert seen["optical"] == seen["thermal"] == {(20.0, 0.003)}
    taller, seen = _urban_case(tmp_path, monkeypatch, "0.003", height="40")
    assert seen["optical"] == seen["thermal"] == {(40.0, 0.003)}
    assert taller != results["0.003"][0]


def test_missing_lai_uses_the_thermal_default(tmp_path):
    lines = (ROOT / "input.csv").read_text(encoding="utf-8-sig").splitlines()
    (tmp_path / "no_lai.csv").write_text(
        "\n".join(line for line in lines if not line.startswith("vegetation_structure,lai,")), encoding="utf-8")
    reference = _rows(run.run(ROOT / "input.csv", tmp_path / "ref" / "output.csv", vza="0;30", vaa="0"))
    missing = _rows(run.run(tmp_path / "no_lai.csv", tmp_path / "missing" / "output.csv", vza="0;30", vaa="0"))
    assert [(row["reflectance_brf"], row["radiance"]) for row in missing] == \
        [(row["reflectance_brf"], row["radiance"]) for row in reference]
    settings = dict(wavelengths=[670.0], leaf_reflectance=0.05, leaf_transmittance=0.03,
                    soil_reflectance=0.12, surface_model="plane", vegetation_model="hom")
    with pytest.raises(ValueError, match="'lai'"):
        reflectance(settings, 670.0, 0.0, 0.0, 30.0, 0.0)


def test_default_input_has_no_warnings(tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        read_input(ROOT / "input.csv")
    _, messages = _run_quiet(ROOT / "input.csv", tmp_path / "output.csv")
    assert messages == []


def test_unknown_and_short_rows_warn_with_line_numbers(tmp_path):
    path = tmp_path / "typo.csv"
    path.write_text("# header\nvegetation_structure,lia,2,typo\ngeometry,vza\n", encoding="utf-8")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        config = read_input(path)
    messages = [str(item.message) for item in caught]
    assert len(messages) == 2
    assert ":2:" in messages[0] and "'lia'" in messages[0]
    assert ":3:" in messages[1]
    assert config["lia"] == "2" and "vza" not in config


def test_non_utf8_input_is_reported(tmp_path):
    path = tmp_path / "gbk.csv"
    path.write_bytes("surface,surface_model,plane,平面地表\n".encode("gbk"))
    with pytest.raises(ValueError, match="UTF-8"):
        read_input(path)


def test_numeric_errors_name_the_key_and_value(tmp_path):
    with pytest.raises(ValueError, match=r"lai.*'abc'"):
        run.run(ROOT / "input.csv", tmp_path / "output.csv", lai="abc")
    with pytest.raises(ValueError, match=r"soil_sunlit_c.*'hot'"):
        run.run(ROOT / "input.csv", tmp_path / "output.csv", soil_sunlit_c="hot", wavelengths="10.5",
                leaf_reflectance="0.015", leaf_transmittance="0", soil_reflectance="0.045",
                street_reflectance="0.045", wall_reflectance="0.07", roof_reflectance="0.05")


def test_scenario_rejects_separate_terrain_emissivity():
    from base.scenario_runner import validate_case
    case = dict(name="t", domain="thermal", surface_model="terrain", vegetation_model="bare",
                wavelengths=[10.5], lai_values=[0], signed_vza=[0], angle_wavelength=10.5,
                sza=30.0, saa=0.0, spectrum_vza=0.0, spectrum_vaa=0.0,
                materials={"emissivity_soil": [0.95], "emissivity_terrain": [0.90]},
                temperatures={"soil_sunlit_c": 45.0})
    with pytest.raises(ValueError, match="emissivity_terrain"):
        validate_case(case)
    case["materials"]["emissivity_terrain"] = [0.95]
    validate_case(case)
