"""Local GUI server for the TiRT directional simulator."""

from __future__ import annotations

import csv
import io
import json
import shutil
import sys
import tempfile
import threading
import zipfile
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
from base.inputio import read_input  # noqa: E402
from base.version import __version__  # noqa: E402
from gui.tirteb_server import _check_paths, _time_request, _json_safe  # noqa: E402


def _prepare_overrides(payload: dict, work_dir: Path) -> dict:
    raw_overrides = payload.get("overrides", {})
    if not isinstance(raw_overrides, dict):
        raise ValueError("overrides must be an object")
    overrides = {str(key): value for key, value in raw_overrides.items()}
    geometry_text = payload.get("geometry_text")
    if geometry_text:
        if not isinstance(geometry_text, str):
            raise ValueError("geometry_text must contain text")
        mode = int(float(overrides.get("geometry_mode", 2)))
        source = str(overrides.get("additional_geometry_source", "none")).strip().lower()
        # Ignore retained uploads when the current selection is manual or principal-only.
        if mode == 2 or (mode == 3 and source in {"file", "txt", "read"}):
            geometry_file = work_dir / "geometry.txt"
            geometry_file.write_text(geometry_text, encoding="utf-8")
            if mode == 3:
                overrides["additional_geometry_file"] = str(geometry_file)
            else:
                overrides["geometry_file"] = str(geometry_file)
                overrides["geometry_mode"] = 2
    _check_paths(overrides, work_dir)
    return overrides


def _export_inputs(overrides: dict, work_dir: Path) -> bytes:
    from base.inputio import write_input_snapshot

    input_path = PROJECT_ROOT / "input.csv"
    config = read_input(input_path)
    config.update({key: str(item) for key, item in overrides.items() if item is not None})
    export_dir = work_dir / "export"
    write_input_snapshot(config, input_path, export_dir / "input.csv")
    (export_dir / "README.txt").write_text(
        "Extract this ZIP, then run the existing TiRT runner with this input:\n"
        "python <TiRT-project>/run.py --input <extracted-folder>/input.csv\n"
        "Keep input.csv and inputs/ together. The bundled files use relative paths.\n",
        encoding="utf-8",
    )
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for file_path in sorted(export_dir.rglob("*")):
            if file_path.is_file():
                archive.write(file_path, file_path.relative_to(export_dir).as_posix())
    return buffer.getvalue()


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

    def log_message(self, format, *args):
        """Windowed EXE builds do not have a usable stderr stream."""
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

    def _trusted_request(self) -> bool:
        """Block DNS rebinding and cross-site requests to the local API."""
        bound, port = self.server.server_address[:2]
        names = {"127.0.0.1", "localhost", "[::1]", f"[{bound}]" if ":" in bound else bound}
        hosts = {f"{name}:{port}".lower() for name in names}
        if self.headers.get("Host", "").strip().lower() not in hosts:
            return False
        origin = self.headers.get("Origin")
        return origin is None or origin.strip().lower() in {f"http://{host}" for host in hosts}

    def do_GET(self):  # noqa: N802
        if urlparse(self.path).path.startswith("/api/") and not self._trusted_request():
            self._json(403, {"ok": False, "error": "Forbidden: request must come from the local TiRT page"})
            return
        if urlparse(self.path).path == "/api/version":
            info_path = PROJECT_ROOT / 'build_info.json'
            info = json.loads(info_path.read_text(encoding='utf-8')) if info_path.exists() else {}
            self._json(200, dict(info, version=__version__, frozen=bool(getattr(sys, 'frozen', False))))
            return
        if urlparse(self.path).path == "/api/health":
            self._json(200, {"ok": True, "project": str(PROJECT_ROOT), "engine": "tirt", "version": __version__})
            return
        super().do_GET()

    def do_POST(self):  # noqa: N802
        path = urlparse(self.path).path
        if path.startswith("/api/") and not self._trusted_request():
            self._json(403, {"ok": False, "error": "Forbidden: request must come from the local TiRT page"})
            return
        if path == "/api/shutdown":
            self._json(200, {"ok": True, "message": "TiRT is shutting down"})
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        if path not in {"/api/run", "/api/export", "/api/run_time"}:
            self._json(404, {"error": "Not found"})
            return
        work_dir = None
        try:
            length = int(self.headers.get("Content-Length", "0"))
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            if not isinstance(payload, dict):
                raise ValueError("Request body must be an object")
            work_dir = Path(tempfile.mkdtemp(prefix="tirt-gui-"))
            if path == "/api/run_time":
                self._json(200, _time_request(payload, work_dir))
                return
            overrides = _prepare_overrides(payload, work_dir)
            if path == "/api/export":
                body = _export_inputs(overrides, work_dir)
                self.send_response(200)
                self.send_header("Content-Type", "application/zip")
                self.send_header("Content-Disposition", 'attachment; filename="tirt_input.zip"')
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Cache-Control", "no-store")
                self.end_headers()
                self.wfile.write(body)
                return
            lai_text = str(overrides.pop("lai", "")).strip()
            lai_values = [float(item) for item in lai_text.replace(",", ";").split(";") if item.strip()]
            if not lai_values:
                raise ValueError("LAI must contain at least one positive number")
            if any(not 0 <= item <= 20 for item in lai_values):
                raise ValueError("LAI values must be in 0..20")
            rows = []
            omitted = 0
            for index, lai in enumerate(lai_values):
                report = {}
                output = work_dir / f"output_lai_{index}.csv"
                run(PROJECT_ROOT / "input.csv", output, report=report,
                    **dict(overrides, lai=str(lai)))
                omitted += int(report.get("omitted_backface", 0))
                for row in _read_csv(output):
                    row["lai"] = lai
                    rows.append(row)
            self._json(200, {"ok": True, "rows": rows, "count": len(rows),
                             "lai_count": len(lai_values), "omitted_backface": omitted})
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
