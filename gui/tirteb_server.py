"""Local GUI server for the TIRTEB temporal and directional simulators."""

from __future__ import annotations

import csv
import json
import math
import shutil
import sys
import tempfile
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse


_BUNDLE_ROOT = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
GUI_ROOT = _BUNDLE_ROOT / "gui" if (_BUNDLE_ROOT / "gui").exists() else Path(__file__).resolve().parent
PROJECT_ROOT = GUI_ROOT.parent
sys.path.insert(0, str(PROJECT_ROOT))

from angle_run import run_angle  # noqa: E402
from tirteb_run import run  # noqa: E402

TIRTEB_INPUT = PROJECT_ROOT / "input_tirteb.csv"


_PATH_KEYS = {
    "meteo", "esun_file", "esky_file", "optical_coefficients_file",
    "directions_file", "spectrals_file", "vegetation_spectrum_file",
    "soil_spectrum_file", "urban_spectrum_file", "urban_roof_spectrum_file",
    "urban_wall_spectrum_file", "urban_street_spectrum_file", "geometry_file",
    "terrain_spectrum_file", "roof_spectrum_file", "wall_spectrum_file",
    "street_spectrum_file", "additional_geometry_file", "observation_geometry_file",
}

def _json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value



def _write_run_input(source: Path, target: Path, overrides: dict[str, object]) -> None:
    with source.open(newline="", encoding="utf-8-sig") as file:
        rows = list(csv.reader(file))
    values = {str(key).lower(): value for key, value in overrides.items() if not str(key).startswith("_")}
    for row in rows:
        if len(row) < 3:
            continue
        key = row[1].strip().lower()
        changed = key in values
        value = values.pop(key) if changed else row[2]
        if key in _PATH_KEYS and value not in (None, ""):
            file_path = Path(str(value))
            if not file_path.is_absolute():
                value = str((PROJECT_ROOT / file_path).resolve())
        if changed or key in _PATH_KEYS:
            row[2] = str(value)
    for key, value in values.items():
        if value not in (None, ""):
            if key in _PATH_KEYS:
                path = Path(str(value))
                if not path.is_absolute():
                    value = str((PROJECT_ROOT / path).resolve())
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


def _angle_request(payload: dict, work_dir: Path) -> dict:
    """Build a temporary angle input and run the directional simulator."""
    overrides = {str(key): value for key, value in payload.get("overrides", {}).items()}
    overrides["spectral_source"] = "direct"

    geometry_text = payload.get("geometry_text")
    if geometry_text:
        geometry_path = work_dir / "geometry.txt"
        geometry_path.write_text(str(geometry_text), encoding="utf-8")
        overrides["additional_geometry_file"] = str(geometry_path)
        overrides["geometry_file"] = str(geometry_path)
        overrides["additional_geometry_source"] = "file"

    uploaded_files = payload.get("files", {}) or {}
    file_overrides = {
        "wavelengths": ("spectral_source", "file", "spectrals_file", "wavelengths.txt"),
        "leaf_spectrum": ("vegetation_spectrum_source", "file", "vegetation_spectrum_file", "leaf_spectrum.txt"),
        "soil_spectrum": ("soil_spectrum_source", "file", "soil_spectrum_file", "soil_spectrum.txt"),
        "terrain_spectrum": ("terrain_spectrum_source", "file", "terrain_spectrum_file", "terrain_spectrum.txt"),
        "roof_spectrum": ("roof_spectrum_source", "file", "roof_spectrum_file", "roof_spectrum.txt"),
        "wall_spectrum": ("wall_spectrum_source", "file", "wall_spectrum_file", "wall_spectrum.txt"),
        "street_spectrum": ("street_spectrum_source", "file", "street_spectrum_file", "street_spectrum.txt"),
    }
    for key, (source_key, source_value, path_key, filename) in file_overrides.items():
        contents = uploaded_files.get(key)
        if not contents:
            continue
        file_path = work_dir / filename
        file_path.write_text(str(contents), encoding="utf-8")
        overrides[source_key] = source_value
        overrides[path_key] = str(file_path)

    input_path = work_dir / "input.csv"
    _write_run_input(TIRTEB_INPUT, input_path, overrides)
    output_path = work_dir / "output.csv"
    run_angle(input_path, output_path)
    rows = _read_csv(output_path)
    return {"ok": True, "rows": rows, "count": len(rows)}


def _time_overrides(payload: dict) -> dict:
    """Translate GUI geometry/material names to the temporal solver contract."""
    if not isinstance(payload, dict):
        raise ValueError("Request body must be an object")
    raw, files = payload.get("overrides", {}), payload.get("files") or {}
    if not isinstance(raw, dict) or not isinstance(files, dict):
        raise ValueError("overrides and files must be objects")
    overrides = {str(key): value for key, value in raw.items()}
    for gui_key, solver_key in {
        "canopy_height": "height", "slope_angle": "terrain_slope",
        "slope_aspect": "terrain_aspect", "crown_rad_a": "crown_width",
        "emissivity_roof": "urban_emissivity_roof", "emissivity_wall": "urban_emissivity_wall",
        "emissivity_street": "urban_emissivity_street",
    }.items():
        if gui_key in overrides:
            overrides[solver_key] = overrides[gui_key]
    # Crown height in TIRTEB is the ellipsoid vertical semi-axis, not trunk height.
    if "crown_rad_b" in overrides:
        overrides["crown_height"] = overrides["crown_rad_b"]
    if "crown_rad_a" in overrides:
        overrides["crown_radius"] = overrides["crown_rad_a"]
    for component in ("leaf", "soil"):
        key = f"emissivity_{component}"
        if key in overrides:
            values = [float(item) for item in str(overrides[key]).replace(",", ";").split(";")]
            if not values or any(not math.isfinite(e) or not 0 <= e <= 1 for e in values):
                raise ValueError(f"{key} must be between 0 and 1")
            overrides[f"{component}_reflectance"] = ";".join(str(1 - e) for e in values)
            if component == "leaf":
                overrides["leaf_transmittance"] = "0"
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
        "wavelengths": ("observation_spectral_mode", "file", "spectrals_file", "observation_spectrals.txt"),
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
        if key == "wavelengths":
            overrides["spectral_source"] = "file"
            overrides["observation_spectral_mode"] = "file"

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
        "series": series,
        "observations": observations,
        "intermediate": intermediate,
        "count": len(series),
        "observation_count": len(observations),
        "intermediate_count": len(intermediate),
    }


class Handler(SimpleHTTPRequestHandler):
    """Serve both TIRTEB UIs and their local simulation endpoints."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(GUI_ROOT), **kwargs)

    def log_message(self, format, *args):
        """Keep windowed builds independent of unavailable stdout/stderr."""
        return

    def _json(self, status: int, payload: dict) -> None:
        body = json.dumps(_json_safe(payload), ensure_ascii=False, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def end_headers(self):
        if not self.path.startswith("/api/"):
            self.send_header("Cache-Control", "no-store, no-cache, must-revalidate")
        super().end_headers()

    def do_GET(self):  # noqa: N802
        path = urlparse(self.path).path
        if path == "/api/health":
            self._json(200, {"ok": True, "project": str(PROJECT_ROOT), "engine": "tirteb"})
            return
        if path == "/":
            self.send_response(302)
            self.send_header("Location", "/time.html")
            self.end_headers()
            return
        super().do_GET()

    def do_POST(self):  # noqa: N802
        path = urlparse(self.path).path
        if path not in {"/api/run_time", "/api/run_angle", "/api/shutdown"}:
            self._json(404, {"error": "Not found"})
            return

        if path == "/api/shutdown":
            self._json(200, {"ok": True})
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return

        work_dir = None
        try:
            length = int(self.headers.get("Content-Length", "0"))
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            overrides = {str(key): value for key, value in payload.get("overrides", {}).items()}
            work_dir = Path(tempfile.mkdtemp(prefix="tirteb-gui-"))

            if path == "/api/run_angle":
                self._json(200, _angle_request(payload, work_dir))
                return

            self._json(200, _time_request(payload, work_dir))
        except Exception as exc:
            self._json(400, {"ok": False, "error": str(exc)})
        finally:
            if work_dir is not None:
                shutil.rmtree(work_dir, ignore_errors=True)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Run the TIRTEB temporal GUI")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8901)
    args = parser.parse_args()
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"TIRTEB GUI: http://{args.host}:{args.port}/time.html")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping TIRTEB GUI")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
