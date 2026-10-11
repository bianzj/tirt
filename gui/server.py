"""Local GUI server for all TiRT pages (directional, energy balance and multi-angle)."""

from __future__ import annotations

import shutil
import sys
import tempfile
from http.server import ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gui.tirteb_server import (  # noqa: E402
    ANGLE_UPLOADS,
    PROJECT_ROOT,
    Handler as EnergyBalanceHandler,
    _read_csv,
    apply_uploads,
    client_overrides,
)
from run import run  # noqa: E402


def _directional_request(payload: dict, work_dir: Path) -> dict:
    overrides = client_overrides(payload)
    apply_uploads(payload.get("files") or {}, ANGLE_UPLOADS, work_dir, overrides)
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
    return {"ok": True, "rows": rows, "count": len(rows)}


class Handler(EnergyBalanceHandler):
    """Serve all TiRT pages and simulation endpoints from one local server."""

    def do_GET(self):  # noqa: N802
        if urlparse(self.path).path == "/api/health":
            if self.request_rejected():
                return
            self._json(200, {
                "ok": True,
                "engine": "tirt",
                "features": ["directional", "energy_balance", "multi_angle"],
            })
            return
        super().do_GET()

    def do_POST(self):  # noqa: N802
        if urlparse(self.path).path != "/api/run":
            super().do_POST()
            return
        if self.request_rejected():
            return
        work_dir = Path(tempfile.mkdtemp(prefix="tirt-gui-"))
        try:
            self._json(200, _directional_request(self.read_json_body(), work_dir))
        except Exception as exc:
            self._json(400, {"ok": False, "error": str(exc)})
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Run the TiRT GUI server")
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
