"""Run the TiRT local HTTP service.

The desktop window is owned by Electron.  This process only owns the local
simulation service, which makes shutdown and process-tree cleanup explicit.
"""

from __future__ import annotations

import argparse
import ctypes
import errno
import json
import sys
from http.server import ThreadingHTTPServer
from pathlib import Path

try:
    from .server import Handler
except ImportError:
    from server import Handler


PORT_SEARCH_COUNT = 20


def _show_error(message: str) -> None:
    """Report startup failures even when running as a windowed executable."""
    if sys.platform == "win32":
        ctypes.windll.user32.MessageBoxW(None, message, "TiRT", 0x10)
    else:
        print(message, file=sys.stderr)


def _create_server(host: str, first_port: int) -> tuple[ThreadingHTTPServer, int]:
    last_error: OSError | None = None
    for port in range(first_port, first_port + PORT_SEARCH_COUNT):
        try:
            server = ThreadingHTTPServer((host, port), Handler)
            return server, int(server.server_address[1])
        except OSError as exc:
            last_error = exc
            if exc.errno not in {errno.EADDRINUSE, 10048}:
                raise
    last_port = first_port + PORT_SEARCH_COUNT - 1
    raise OSError(f"Ports {first_port}-{last_port} are unavailable") from last_error


def _write_ready_file(path: str | None, host: str, port: int) -> None:
    if not path:
        return
    ready_path = Path(path)
    ready_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = ready_path.with_suffix(f"{ready_path.suffix}.tmp")
    temporary.write_text(json.dumps({"ok": True, "host": host, "port": port}), encoding="utf-8")
    temporary.replace(ready_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the TiRT local simulation service")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-browser", action="store_true", help="兼容旧启动脚本：只启动服务")
    parser.add_argument("--server-only", action="store_true", help="只启动服务，由桌面壳层管理窗口")
    parser.add_argument("--ready-file", help="服务就绪后写入 JSON 文件")
    args = parser.parse_args()

    server = None
    try:
        server, port = _create_server(args.host, args.port)
        server.daemon_threads = True
        url = f"http://{args.host}:{port}/"
        _write_ready_file(args.ready_file, args.host, port)
        if sys.stdout is not None:
            print(f"TiRT service: {url}", flush=True)
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    except Exception as exc:
        _show_error(f"TiRT service could not start.\n\n{exc}")
    finally:
        if server is not None:
            server.server_close()


if __name__ == "__main__":
    main()
