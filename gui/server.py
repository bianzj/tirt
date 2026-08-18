"""Local GUI server for the TiRT Three.js simulator."""

from __future__ import annotations

import csv
import json
import os
import shutil
import sys
import tempfile
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse


<<<<<<< Updated upstream
def _project_root() -> Path:
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS).resolve()
    return Path(__file__).resolve().parent.parent


PROJECT_ROOT = _project_root()
GUI_ROOT = PROJECT_ROOT / "gui"
=======
GUI_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = GUI_ROOT.parent
TIRTEB_ROOT = Path(os.environ.get("TIRTEB_ROOT", str(PROJECT_ROOT.parent / "tirteb"))).resolve()
>>>>>>> Stashed changes
sys.path.insert(0, str(PROJECT_ROOT))
os.environ.setdefault("MPLBACKEND", "Agg")

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
    """Serve the UI and expose a small JSON endpoint for one simulation."""

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
            self._json(200, {"ok": True, "project": str(PROJECT_ROOT)})
            return
        super().do_GET()

    def do_POST(self):  # noqa: N802
        path = urlparse(self.path).path
        if path not in {"/api/run", "/api/run_time"}:
            self._json(404, {"error": "Not found"})
            return
        if path == "/api/run_time":
            self._run_time_series()
            return
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
                if contents:
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
            shutil.rmtree(work_dir, ignore_errors=True)
        except Exception as exc:  # Keep the browser error actionable.
            self._json(400, {"ok": False, "error": str(exc)})

    def _run_time_series(self) -> None:
        """Run the TIRTEB time-step engine in an isolated Python process.

        TIRTEB and TiRT both have a ``base`` package, so importing both engines
        into this server would make Python resolve modules from the wrong
        project. A subprocess keeps each engine's package namespace isolated.
        """
        work_dir = None
        try:
            if not (TIRTEB_ROOT / "run.py").is_file() or not (TIRTEB_ROOT / "input.csv").is_file():
                raise FileNotFoundError(
                    f"未找到 tirteb 项目: {TIRTEB_ROOT}。可设置环境变量 TIRTEB_ROOT 指向 tirteb 根目录。"
                )
            length = int(self.headers.get("Content-Length", "0"))
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            overrides = {str(key): value for key, value in payload.get("overrides", {}).items()}
            work_dir = Path(tempfile.mkdtemp(prefix="tirteb-gui-"))
            meteo_text = payload.get("meteo_text")
            if meteo_text:
                meteo_path = work_dir / "meteo.txt"
                meteo_path.write_text(str(meteo_text), encoding="utf-8")
                overrides["meteo"] = str(meteo_path)

            uploaded_files = payload.get("files", {}) or {}
            time_file_overrides = {
                "directions": ("observation_direction_mode", "file", "directions_file", "observation_directions.txt"),
                "wavelengths": ("observation_spectral_mode", "file", "spectrals_file", "observation_spectrals.txt"),
                "vegetation_spectrum": ("vegetation_spectrum_source", "measured", "vegetation_spectrum_file", "leaf_spectrum.txt"),
                "soil_spectrum": ("soil_spectrum_source", "measured", "soil_spectrum_file", "soil_spectrum.txt"),
                "urban_spectrum": ("urban_spectrum_source", "measured", "urban_spectrum_file", "urban_spectrum.txt"),
            }
            for key, (source_key, source_value, path_key, filename) in time_file_overrides.items():
                contents = uploaded_files.get(key)
                if contents:
                    file_path = work_dir / filename
                    file_path.write_text(str(contents), encoding="utf-8")
                    overrides[source_key] = source_value
                    overrides[path_key] = str(file_path)
                    if key == "directions":
                        overrides["direction_source"] = "file"
                    if key == "wavelengths":
                        overrides["spectral_source"] = "file"
                        overrides["observation_spectral_mode"] = "file"

            output_path = work_dir / "output.csv"
            call_file = work_dir / "call.json"
            call_file.write_text(json.dumps({
                "input": str(TIRTEB_ROOT / "input.csv"),
                "output": str(output_path),
                "overrides": overrides,
            }, ensure_ascii=False), encoding="utf-8")
            runner = (
                "import json, sys; "
                "from pathlib import Path; "
                "from run import run; "
                "p=json.loads(Path(sys.argv[1]).read_text(encoding='utf-8')); "
                "run(p['input'], p['output'], **p['overrides'])"
            )
            import subprocess
            completed = subprocess.run(
                [sys.executable, "-c", runner, str(call_file)],
                cwd=str(TIRTEB_ROOT),
                capture_output=True,
                text=True,
                timeout=300,
                check=False,
            )
            if completed.returncode != 0:
                detail = (completed.stderr or completed.stdout).strip()
                raise RuntimeError(detail or "tirteb 时间序列运行失败")

            series = _read_csv(output_path)
            observations_path = output_path.with_name("observations.csv")
            observations = _read_csv(observations_path) if observations_path.exists() else []
            self._json(200, {
                "ok": True,
                "series": series,
                "observations": observations,
                "count": len(series),
                "observation_count": len(observations),
            })
        except Exception as exc:  # Keep the browser error actionable.
            self._json(400, {"ok": False, "error": str(exc)})
        finally:
            if work_dir is not None:
                shutil.rmtree(work_dir, ignore_errors=True)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Run the TiRT Three.js GUI")
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
