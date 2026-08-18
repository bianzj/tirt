"""Local GUI server for the TiRT directional simulator."""

from __future__ import annotations

import csv
import json
import shutil
import sys
import tempfile
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse


def _project_root() -> Path:
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS).resolve()
    return Path(__file__).resolve().parent.parent


PROJECT_ROOT = _project_root()
GUI_ROOT = PROJECT_ROOT / "gui"
sys.path.insert(0, str(PROJECT_ROOT))

from run import run  # noqa: E402


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


class Handler(SimpleHTTPRequestHandler):
    """Serve the TiRT directional UI and its local simulation endpoint."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(GUI_ROOT), **kwargs)

    def _json(self, status: int, payload: dict) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):  # noqa: N802
        if urlparse(self.path).path == "/api/health":
            self._json(200, {"ok": True, "project": str(PROJECT_ROOT), "engine": "tirt"})
            return
        super().do_GET()

    def do_POST(self):  # noqa: N802
        if urlparse(self.path).path != "/api/run":
            self._json(404, {"error": "Not found"})
            return
        work_dir = None
        try:
            length = int(self.headers.get("Content-Length", "0"))
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            overrides = {str(key): value for key, value in payload.get("overrides", {}).items()}
            work_dir = Path(tempfile.mkdtemp(prefix="tirt-gui-"))
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
                if key == "wavelengths":
                    overrides["spectral_source"] = "file"
                    overrides["observation_spectral_mode"] = "file"
            geometry_text = payload.get("geometry_text")
            if geometry_text:
                geometry_file = work_dir / "geometry.txt"
                geometry_file.write_text(str(geometry_text), encoding="utf-8")
                if str(overrides.get("geometry_mode", "")) == "3":
                    overrides["additional_geometry_file"] = str(geometry_file)
                    overrides["additional_geometry_source"] = "file"
                else:
                    overrides["geometry_file"] = str(geometry_file)
                    overrides["geometry_mode"] = 2
            output = run(PROJECT_ROOT / "input.csv", work_dir / "output.csv", **overrides)
            rows = _read_csv(output)
            self._json(200, {"ok": True, "rows": rows, "count": len(rows)})
        except Exception as exc:
            self._json(400, {"ok": False, "error": str(exc)})
        finally:
            if work_dir is not None:
                shutil.rmtree(work_dir, ignore_errors=True)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Run the TiRT directional GUI")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"TiRT GUI: http://{args.host}:{args.port}/")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping TiRT GUI")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
