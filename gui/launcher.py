"""Executable launcher for the STREAM TIRT local GUI."""

from __future__ import annotations

import argparse
import ctypes
import errno
import json
import sys
import threading
import traceback
import webbrowser
from pathlib import Path
from http.server import ThreadingHTTPServer

try:
    from .server import Handler, PROJECT_ROOT
except ImportError:
    from server import Handler, PROJECT_ROOT
from base.version import __version__


PORT_SEARCH_COUNT = 20


def _show_error(message: str) -> None:
    """Report startup failures even when running as a windowed executable."""
    if sys.platform == "win32":
        ctypes.windll.user32.MessageBoxW(None, message, "STREAM TIRT", 0x10)
    else:
        print(message, file=sys.stderr)


def _create_server(host: str, first_port: int) -> tuple[ThreadingHTTPServer, int]:
    last_error: OSError | None = None
    for port in range(first_port, first_port + PORT_SEARCH_COUNT):
        try:
            return ThreadingHTTPServer((host, port), Handler), port
        except OSError as exc:
            last_error = exc
            if exc.errno not in {errno.EADDRINUSE, 10048}:
                raise
    last_port = first_port + PORT_SEARCH_COUNT - 1
    raise OSError(f"Ports {first_port}-{last_port} are unavailable") from last_error


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the STREAM TIRT GUI")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument("--version", action="version", version=f"TiRT {__version__}")
    parser.add_argument("--verify-model", type=Path,
                        help="Recompute bundled model validation and save results to this directory.")
    args = parser.parse_args()

    if args.verify_model:
        # This runs the same numerical replay used for the published source.
        output = args.verify_model.resolve()
        output.mkdir(parents=True, exist_ok=True)
        try:
            from scripts.validate_current_model import main as verify_model
            sys.argv = [sys.argv[0], '--bundle', str(PROJECT_ROOT / 'docs/validation/current_model_20261021'),
                        '--output', str(output)]
            verify_model()
        except Exception:
            (output / 'error.json').write_text(json.dumps(
                {'ok': False, 'version': __version__, 'traceback': traceback.format_exc()},
                ensure_ascii=False, indent=2), encoding='utf-8')
            raise SystemExit(1)
        return

    server = None
    try:
        server, port = _create_server(args.host, args.port)
        url = f"http://{args.host}:{port}/"
        if args.no_browser:
            if sys.stdout is not None:
                print(f"STREAM TIRT: {url}")
        else:
            opener = threading.Timer(0.5, webbrowser.open, args=(url,), kwargs={"new": 2})
            opener.daemon = True
            opener.start()
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    except Exception as exc:
        _show_error(f"STREAM TIRT could not start.\n\n{exc}")
    finally:
        if server is not None:
            server.server_close()


if __name__ == "__main__":
    main()
