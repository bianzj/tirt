"""Energy-balance (time series) request handling for the TiRT GUI server."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path


_BUNDLE_ROOT = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
GUI_ROOT = _BUNDLE_ROOT / "gui" if (_BUNDLE_ROOT / "gui").exists() else Path(__file__).resolve().parent
PROJECT_ROOT = GUI_ROOT.parent
sys.path.insert(0, str(PROJECT_ROOT))

from tirteb_run import INPUT_FILE_KEYS, run  # noqa: E402

TIRTEB_INPUT = PROJECT_ROOT / "input_tirteb.csv"


_PATH_KEYS = {
    "meteo", "esun_file", "esky_file", "optical_coefficients_file",
    "directions_file", "spectrals_file", "vegetation_spectrum_file",
    "soil_spectrum_file", "urban_spectrum_file", "urban_roof_spectrum_file",
    "urban_wall_spectrum_file", "urban_street_spectrum_file", "geometry_file",
    "terrain_spectrum_file", "roof_spectrum_file", "wall_spectrum_file",
    "street_spectrum_file", "additional_geometry_file", "observation_geometry_file",
} | set(INPUT_FILE_KEYS)


def _normal_key(key) -> str:
    """Canonical parameter name, identical to what base.inputio.read_input sees."""
    return str(key).strip().lower()


def _is_path_key(name: str) -> bool:
    return name in _PATH_KEYS or name.endswith(("_file", "_dir", "_path"))


def _absolute(value) -> str:
    path = Path(str(value))
    return str((path if path.is_absolute() else PROJECT_ROOT / path).resolve())

def _json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _check_paths(overrides: dict, work_dir: Path) -> None:
    """Client-supplied paths may only name bundled data or this request's uploads."""
    roots = ((PROJECT_ROOT / "data").resolve(), Path(work_dir).resolve())
    seen = set()
    for key, value in overrides.items():
        name = _normal_key(key)
        if not name or name.startswith("_"):
            continue
        if name in seen:
            raise ValueError(f"Duplicate parameter after normalisation: {key!r}")
        seen.add(name)
        if any(ch in name for ch in ",\r\n\"") or any(ch in str(value) for ch in "\r\n"):
            raise ValueError(f"Invalid characters in parameter {key!r}")
        if value in (None, "") or not _is_path_key(name):
            continue
        path = Path(_absolute(value))
        if not any(path.is_relative_to(root) for root in roots):
            raise ValueError(f"{key} must be a file inside the project data folder: {value}")


def _write_run_input(source: Path, target: Path, overrides: dict[str, object]) -> None:
    with source.open(newline="", encoding="utf-8-sig") as file:
        rows = list(csv.reader(file))
    # Normalise exactly like read_input (strip + lower) so the checked key is the used key.
    values = {_normal_key(key): value for key, value in overrides.items()
              if _normal_key(key) and not _normal_key(key).startswith("_")}
    for row in rows:
        if len(row) < 3:
            continue
        key = _normal_key(row[1])
        changed = key in values
        value = values.pop(key) if changed else row[2]
        if _is_path_key(key) and str(value).strip() not in ("",):
            value = _absolute(str(value).strip())
        if changed or _is_path_key(key):
            row[2] = str(value)
    for key, value in values.items():
        if value not in (None, ""):
            if _is_path_key(key):
                value = _absolute(value)
            rows.append(["gui", key, str(value), "GUI run value"])
    with target.open("w", newline="", encoding="utf-8") as file:
        csv.writer(file, lineterminator="\n").writerows(rows)


def _read_csv(path: Path) -> list[dict[str, float | str]]:
    with path.open(newline="", encoding="utf-8-sig") as file:
        rows = []
        for row in csv.DictReader(file):
            converted: dict[str, float | str] = {}
            for key, value in row.items():
                if value in (None, ""):
                    continue
                try:
                    converted[key] = float(value)
                except ValueError:
                    converted[key] = value
            rows.append(converted)
        return rows


def _time_overrides(payload: dict) -> dict:
    """Translate GUI geometry/material names to the temporal solver contract."""
    if not isinstance(payload, dict):
        raise ValueError("Request body must be an object")
    raw, files = payload.get("overrides", {}), payload.get("files") or {}
    if not isinstance(raw, dict) or not isinstance(files, dict):
        raise ValueError("overrides and files must be objects")
    overrides = {}
    for key, value in raw.items():
        name = _normal_key(key)
        if name in overrides:
            raise ValueError(f"Duplicate parameter after normalisation: {key!r}")
        overrides[name] = value
    for gui_key, solver_key in {
        "canopy_height": "height", "slope_angle": "terrain_slope",
        "slope_aspect": "terrain_aspect",
    }.items():
        if gui_key in overrides:
            overrides[solver_key] = overrides.pop(gui_key)
    if overrides.get("vegetation_model") == "row" and "row_height" in overrides:
        overrides["height"] = overrides["row_height"]
    if any(key in overrides for key in ("vza", "vaa", "geometry_mode")):
        overrides["direction_source"] = "direct"
    return overrides


def _time_request(payload: dict, work_dir: Path) -> dict:
    overrides = _time_overrides(payload)
    meteo_text = payload.get("meteo_text")
    if meteo_text:
        meteo_path = work_dir / "meteo.txt"
        meteo_path.write_text(str(meteo_text), encoding="utf-8")
        overrides["meteo"] = str(meteo_path)

    uploaded_files = dict(payload.get("files", {}) or {})
    # Accept cached clients using the earlier upload names.
    for old, current in (("leaf_optics", "vegetation_spectrum"), ("soil_optics", "soil_spectrum")):
        if uploaded_files.get(old) and not uploaded_files.get(current):
            uploaded_files[current] = uploaded_files[old]
    # The UI supplies pairs; the temporal runner's direct lists are a cross product.
    # Write pairs as a local direction file so no unintended views are generated.
    if not uploaded_files.get("directions") and str(overrides.get("geometry_pair_mode", "")) == "zip":
        def angles(key):
            return [float(x) for x in str(overrides.get(key, "0")).replace(",", ";").split(";")]
        zenith, azimuth = angles("vza"), angles("vaa")
        if len(zenith) != len(azimuth):
            raise ValueError("Each observation zenith must have a matching azimuth")
        if any(not math.isfinite(v) or not 0 <= v < 90 for v in zenith) or any(
                not math.isfinite(v) for v in azimuth):
            raise ValueError("Observation angles require finite azimuths and 0 <= zenith < 90")
        uploaded_files["directions"] = "\n".join(f"{z} {a}" for z, a in zip(zenith, azimuth))
    file_overrides = {
        "directions": ("observation_direction_mode", "file", "directions_file", "observation_directions.txt"),
        "vegetation_spectrum": ("vegetation_spectrum_source", "measured", "vegetation_spectrum_file", "leaf_spectrum.txt"),
        "soil_spectrum": ("soil_spectrum_source", "measured", "soil_spectrum_file", "soil_spectrum.txt"),
        "urban_spectrum": ("urban_spectrum_source", "measured", "urban_spectrum_file", "urban_spectrum.txt"),
    }
    for key, (source_key, source_value, path_key, filename) in file_overrides.items():
        contents = uploaded_files.get(key)
        if not contents:
            continue
        file_path = work_dir / filename
        file_path.write_text(str(contents), encoding="utf-8")
        overrides[source_key] = source_value
        overrides[path_key] = str(file_path)
        if key == "directions":
            overrides["direction_source"] = "file"

    # TiRT-only geometry keys are rejected by the temporal runner.
    for key in [key for key in overrides if _normal_key(key) in {"geometry_mode", "geometry_pair_mode", "geometry_file"}
                or _normal_key(key).startswith(("principal_vza", "additional_geometry_"))]:
        del overrides[key]
    _check_paths(overrides, work_dir)
    input_path = work_dir / "input.csv"
    _write_run_input(TIRTEB_INPUT, input_path, overrides)
    output_path = work_dir / "output.csv"
    output = run(input_path, output_path)
    series = _read_csv(output)
    observations_path = output.with_name("observations.csv")
    observations = _read_csv(observations_path) if observations_path.exists() else []
    intermediate_path = output.with_name("intermediate.csv")
    intermediate = _read_csv(intermediate_path) if intermediate_path.exists() else []
    return {
        "ok": True,
        "scene": {key: overrides[key] for key in ("surface_model", "vegetation_model") if key in overrides},
        "series": series,
        "observations": observations,
        "intermediate": intermediate,
        "count": len(series),
        "observation_count": len(observations),
        "intermediate_count": len(intermediate),
    }
