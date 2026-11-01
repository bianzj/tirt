"""Regression tests for the directional DEM facet-ensemble surface."""
import math
from pathlib import Path

import numpy as np
import pytest

import run
from base import dem


ROOT = Path(__file__).resolve().parents[1]


def _write_ridge_dem(path, n=25, cell=1.0):
    x = np.arange(n) * cell
    # Two triangular ridges running north-south (constant along rows).
    z = np.maximum(0.0, 2.0 - np.abs(x - 8.0) * 0.5) + np.maximum(0.0, 2.0 - np.abs(x - 17.0) * 0.5)
    z = np.tile(z, (n, 1))
    lines = [f"NCOLS {n}", f"NROWS {n}", "XLLCORNER 0", "YLLCORNER 0", f"CELLSIZE {cell:g}", "NODATA_VALUE -9999"]
    lines += [" ".join(f"{value:.5f}" for value in row) for row in z]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _config(dem_path):
    config = run.read_input(ROOT / "input.csv")
    config.update(surface_model="dem", vegetation_model="bare", dem_file=str(dem_path),
                  dem_slope_bin="2", dem_aspect_bin="10", dem_shadow="1",
                  sky_temperature_k="none", vza="0", vaa="0", sza="70", saa="270")
    return config


def _spectra(config):
    values = run._spectra(run._band_table(config, ROOT / "input.csv"))
    return {key: item[2] for key, item in values.items()}


def test_ascii_grid_roundtrip(tmp_path):
    path = tmp_path / "ridges.asc"
    _write_ridge_dem(path)
    config = {"dem_file": str(path)}
    grid = dem.read_grid(config)
    assert grid.z.shape == (25, 25)
    assert grid.cell == pytest.approx(1.0)
    slope, aspect = dem.slope_aspect(grid)
    assert slope.shape == grid.z.shape and np.all(slope >= 0.0)
    assert np.all((aspect >= 0.0) & (aspect < 360.0))


def test_dem_mode_is_canonical_and_radiance_is_finite(tmp_path):
    path = tmp_path / "ridges.asc"
    _write_ridge_dem(path)
    config = _config(path)
    assert run._canonical(config)[0] == "dem"
    config["surface_model"] = "ridge"
    assert run._canonical(config)[0] == "dem"

    spectra = _spectra(config)
    radiance, brightness = run._surface_result(config, "dem", "bare", 10.5, (0.0, 0.0, 70.0, 270.0), spectra)
    assert np.isfinite(radiance) and radiance > 0.0
    assert np.isfinite(brightness) and 250.0 < brightness < 350.0

    plane, plane_bt = run._surface_result(dict(config, surface_model="plane"), "plane", "bare", 10.5,
                                          (0.0, 0.0, 70.0, 270.0), spectra)
    assert brightness < plane_bt, "shaded ridge facets should lower the ensemble brightness temperature"


def test_dem_mode_optical_is_finite(tmp_path):
    path = tmp_path / "ridges.asc"
    _write_ridge_dem(path)
    config = _config(path)
    table = run._band_table(config, ROOT / "input.csv")
    band = {key: float(table[key][0]) for key in run.BAND_PROPERTIES}
    value = run._optical_band(config, "dem", "bare", 0.67, (60.0, 90.0, 70.0, 270.0), band)
    assert np.isfinite(value) and 0.0 < value < 1.0


def test_missing_dem_file_is_a_clear_error(tmp_path):
    config = _config(tmp_path / "missing.asc")
    with pytest.raises(ValueError, match="dem_file"):
        dem.read_grid(config)


def test_csv_entry_resolves_relative_dem_and_writes_output(tmp_path):
    path = tmp_path / "ridges.asc"
    _write_ridge_dem(path)
    template = (ROOT / "input.csv").read_text(encoding="utf-8")
    template = template.replace("surface,surface_model,plane,",
                                "surface,surface_model,dem,")
    template += "\ndem,dem_file,ridges.asc,test grid\ndem,dem_slope_bin,2,test\ndem,dem_aspect_bin,10,test\n"
    input_path = tmp_path / "input.csv"
    input_path.write_text(template, encoding="utf-8")
    output_path = tmp_path / "output.csv"
    run.run(input_path, output_path)
    text = output_path.read_text(encoding="utf-8-sig")
    assert "surface_model" in text and "dem" in text
    assert len(text.splitlines()) > 1
