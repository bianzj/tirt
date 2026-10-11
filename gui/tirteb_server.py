"""Shared request handler for the TiRT local service (energy-balance and multi-angle endpoints).

`gui/server.py` subclasses this handler and is the entry point that serves all pages.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import sys
import tempfile
import threading
from http.server import SimpleHTTPRequestHandler
from pathlib import Path
from urllib.parse import unquote, urlparse


_BUNDLE_ROOT = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
GUI_ROOT = _BUNDLE_ROOT / "gui" if (_BUNDLE_ROOT / "gui").exists() else Path(__file__).resolve().parent
PROJECT_ROOT = GUI_ROOT.parent
sys.path.insert(0, str(PROJECT_ROOT))

from angle_run import run_angle  # noqa: E402
from tirteb_run import run  # noqa: E402

TIRTEB_INPUT = PROJECT_ROOT / "input_tirteb.csv"


_PATH_KEYS = {"meteo"}
_CLIENT_FORBIDDEN_KEYS = {"result_dir", "cases_dir", "case_name", "input_path", "output_path", "model"}
_LOOPBACK_HOSTS = {"127.0.0.1", "localhost", "::1"}
MAX_REQUEST_BYTES = 64 * 1024 * 1024

# upload key -> (source key, source value, path key, file name)
ANGLE_UPLOADS = {
    "wavelengths": ("spectral_source", "file", "spectrals_file", "wavelengths.txt"),
    "leaf_spectrum": ("vegetation_spectrum_source", "file", "vegetation_spectrum_file", "leaf_spectrum.txt"),
    "soil_spectrum": ("soil_spectrum_source", "file", "soil_spectrum_file", "soil_spectrum.txt"),
    "terrain_spectrum": ("terrain_spectrum_source", "file", "terrain_spectrum_file", "terrain_spectrum.txt"),
    "roof_spectrum": ("roof_spectrum_source", "file", "roof_spectrum_file", "roof_spectrum.txt"),
    "wall_spectrum": ("wall_spectrum_source", "file", "wall_spectrum_file", "wall_spectrum.txt"),
    "street_spectrum": ("street_spectrum_source", "file", "street_spectrum_file", "street_spectrum.txt"),
}
TIME_UPLOADS = {
    "directions": ("geometry_mode", "2", "geometry_file", "observation_directions.txt"),
    "wavelengths": ("spectral_source", "file", "spectrals_file", "observation_spectrals.txt"),
    "leaf_optics": ("leaf_optics_source", "file", "leaf_optics_file", "leaf_optics.txt"),
    "soil_optics": ("soil_optics_source", "file", "soil_optics_file", "soil_optics.txt"),
    "urban_spectrum": ("urban_spectrum_source", "file", "urban_spectrum_file", "urban_spectrum.txt"),
}


def _is_path_key(key: str) -> bool:
    return key in _PATH_KEYS or key.endswith("_file")


def _json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def client_overrides(payload: dict) -> dict[str, object]:
    """Return browser-supplied overrides with file paths confined to the project."""
    root = PROJECT_ROOT.resolve()
    overrides: dict[str, object] = {}
    for key, value in (payload.get("overrides") or {}).items():
        key = str(key)
        lowered = key.lower()
        if lowered in _CLIENT_FORBIDDEN_KEYS:
            continue
        if _is_path_key(lowered) and value not in (None, ""):
            path = (root / str(value)).resolve()
            if not path.is_relative_to(root):
                raise ValueError(f"{key} 必须位于项目目录内")
            value = str(path)
        overrides[key] = value
    return overrides


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
        if _is_path_key(key) and value not in (None, ""):
            file_path = Path(str(value))
            if not file_path.is_absolute():
                value = str((PROJECT_ROOT / file_path).resolve())
        if changed or _is_path_key(key):
            row[2] = str(value)
    for key, value in values.items():
        if value not in (None, ""):
            if _is_path_key(key):
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


def apply_uploads(uploaded: dict, table: dict, work_dir: Path, overrides: dict[str, object]) -> list[str]:
    """Write uploaded text files into work_dir and point the overrides at them."""
    applied = []
    for key, (source_key, source_value, path_key, filename) in table.items():
        contents = uploaded.get(key)
        if not contents:
            continue
        file_path = work_dir / filename
        file_path.write_text(str(contents), encoding="utf-8")
        overrides[source_key] = source_value
        overrides[path_key] = str(file_path)
        applied.append(key)
    return applied


def _angle_request(payload: dict, work_dir: Path) -> dict:
    """Build a temporary angle input and run the directional simulator."""
    overrides = client_overrides(payload)
    overrides["spectral_source"] = "direct"

    geometry_text = payload.get("geometry_text")
    if geometry_text:
        geometry_path = work_dir / "geometry.txt"
        geometry_path.write_text(str(geometry_text), encoding="utf-8")
        overrides["additional_geometry_file"] = str(geometry_path)
        overrides["geometry_file"] = str(geometry_path)
        overrides["additional_geometry_source"] = "file"

    apply_uploads(payload.get("files") or {}, ANGLE_UPLOADS, work_dir, overrides)
    input_path = work_dir / "input.csv"
    _write_run_input(TIRTEB_INPUT, input_path, overrides)
    output_path = work_dir / "output.csv"
    run_angle(input_path, output_path)
    rows = _read_csv(output_path)
    return {"ok": True, "rows": rows, "count": len(rows)}


def _time_request(payload: dict, work_dir: Path) -> dict:
    """Build a temporary energy-balance input and run the time series."""
    overrides = client_overrides(payload)
    meteo_text = payload.get("meteo_text")
    if meteo_text:
        meteo_path = work_dir / "meteo.txt"
        meteo_path.write_text(str(meteo_text), encoding="utf-8")
        overrides["meteo"] = str(meteo_path)

    apply_uploads(payload.get("files") or {}, TIME_UPLOADS, work_dir, overrides)

    input_path = work_dir / "input.csv"
    _write_run_input(TIRTEB_INPUT, input_path, overrides)
    output = run(input_path, work_dir / "output.csv")
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
    """Serve the GUI pages plus the energy-balance and multi-angle endpoints."""

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

    def request_rejected(self) -> bool:
        """Reject DNS-rebinding hosts and cross-origin (CSRF) requests to the local service."""
        host = self.headers.get("Host", "")
        bound_host, port = self.server.server_address[:2]
        allowed = True
        if bound_host in _LOOPBACK_HOSTS:
            allowed = host in {f"127.0.0.1:{port}", f"localhost:{port}", f"[::1]:{port}"}
        origin = self.headers.get("Origin")
        if origin is not None and origin != f"http://{host}":
            allowed = False
        if not allowed:
            self._json(403, {"ok": False, "error": "Request rejected"})
        return not allowed

    def read_json_body(self) -> dict:
        length = int(self.headers.get("Content-Length", "0"))
        if length < 0 or length > MAX_REQUEST_BYTES:
            raise ValueError("Invalid request size")
        payload = json.loads(self.rfile.read(length).decode("utf-8")) if length else {}
        if not isinstance(payload, dict):
            raise ValueError("Request body must be a JSON object")
        return payload

    def do_GET(self):  # noqa: N802
        if self.request_rejected():
            return
        if Path(unquote(urlparse(self.path).path)).suffix.lower() in {".py", ".pyc"}:
            self._json(404, {"error": "Not found"})
            return
        super().do_GET()

    def do_POST(self):  # noqa: N802
        if self.request_rejected():
            return
        path = urlparse(self.path).path
        if path == "/api/shutdown":
            self._json(200, {"ok": True})
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        handlers = {"/api/run_time": _time_request, "/api/run_angle": _angle_request}
        if path not in handlers:
            self._json(404, {"error": "Not found"})
            return
        work_dir = Path(tempfile.mkdtemp(prefix="tirt-gui-"))
        try:
            self._json(200, handlers[path](self.read_json_body(), work_dir))
        except Exception as exc:
            self._json(400, {"ok": False, "error": str(exc)})
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)
