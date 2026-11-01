"""Executable launcher for the RTLab TIRT local GUI."""

from __future__ import annotations

import argparse
import ctypes
import errno
import json
import os
import stat
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


def _show_error(message: str, dialog: bool = True) -> None:
    """Report startup failures even when running as a windowed executable.

    The desktop window (``--server-only``) reads stderr and shows its own error
    dialog for a non-zero exit, so no blocking MessageBox is shown there.
    """
    if dialog and sys.platform == "win32":
        ctypes.windll.user32.MessageBoxW(None, message, "RTLab TIRT", 0x10)
        return
    stream = sys.stderr or sys.__stderr__
    if stream is not None:
        try:
            print(message, file=stream, flush=True)
        except (OSError, ValueError):
            pass


def _stdin_is_pipe() -> bool:
    """True when stdin is an anonymous/named pipe (as spawned by the desktop window).

    A console or the null device is not watched: /dev/null or NUL would report EOF
    immediately and stop a manually started ``--server-only`` server.
    """
    try:
        return sys.stdin is not None and stat.S_ISFIFO(os.fstat(sys.stdin.fileno()).st_mode)
    except (OSError, ValueError, AttributeError):
        return False


def _watch_parent(server: ThreadingHTTPServer, stream=None) -> threading.Thread:
    """Stop the server when the parent closes our stdin (normal quit, crash or kill).

    The desktop window spawns the backend with ``stdin: 'pipe'`` and never writes to
    it; the operating system closes the pipe when the parent process ends for any
    reason, so a blocking read returns EOF and the backend does not stay orphaned.
    """
    stream = stream if stream is not None else sys.stdin

    def watch():
        try:
            source = getattr(stream, "buffer", stream)
            while source.read(1):
                pass
        except (OSError, ValueError):
            pass
        server.shutdown()

    thread = threading.Thread(target=watch, name="tirt-parent-watch", daemon=True)
    thread.start()
    return thread


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


def main(argv: list[str] | None = None) -> int:
    """Run the GUI server; return the process exit code (0 = normal stop)."""
    parser = argparse.ArgumentParser(description="Run the RTLab TIRT GUI")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument("--server-only", action="store_true",
                        help="Do not open a browser (used by the desktop window).")
    parser.add_argument("--ready-file", type=Path,
                        help="Write the bound host and port as JSON to this file once the server is listening.")
    parser.add_argument("--version", action="version", version=f"TiRT {__version__}")
    parser.add_argument("--verify-model", type=Path,
                        help="Recompute bundled model validation and save results to this directory.")
    parser.add_argument("--exit-on-stdin-eof", action="store_true",
                        help="Stop when stdin closes (automatic with --server-only when stdin is a pipe).")
    args = parser.parse_args(argv)

    if args.verify_model:
        # This runs the same numerical replay used for the published source.
        output = args.verify_model.resolve()
        output.mkdir(parents=True, exist_ok=True)
        try:
            from scripts.validate_current_model import main as verify_model
            sys.argv = [sys.argv[0], '--bundle', str(PROJECT_ROOT / 'validation/reference'),
                        '--output', str(output)]
            verify_model()
        except Exception:
            (output / 'error.json').write_text(json.dumps(
                {'ok': False, 'version': __version__, 'traceback': traceback.format_exc()},
                ensure_ascii=False, indent=2), encoding='utf-8')
            return 1
        return 0

    server = None
    started = False
    try:
        server, port = _create_server(args.host, args.port)
        port = server.server_address[1]
        url = f"http://{args.host}:{port}/"
        if args.exit_on_stdin_eof or (args.server_only and _stdin_is_pipe()):
            _watch_parent(server)
        if args.ready_file:
            args.ready_file.write_text(json.dumps({"host": args.host, "port": port}), encoding="utf-8")
        if args.no_browser or args.server_only:
            if sys.stdout is not None:
                print(f"RTLab TIRT: {url}", flush=True)
        else:
            opener = threading.Timer(0.5, webbrowser.open, args=(url,), kwargs={"new": 2})
            opener.daemon = True
            opener.start()
        started = True
        server.serve_forever()
    except KeyboardInterrupt:
        return 0
    except Exception as exc:
        stage = "stopped unexpectedly" if started else "could not start"
        _show_error(f"RTLab TIRT {stage}.\n\n{exc}\n\n{traceback.format_exc()}",
                    dialog=not args.server_only)
        return 1
    finally:
        if server is not None:
            server.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
