"""Executable launcher for the TiRT local GUI."""

from __future__ import annotations

import argparse
import threading
import webbrowser
from http.server import ThreadingHTTPServer

from server import Handler


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the TiRT GUI")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    url = f"http://{args.host}:{args.port}/"
    print(f"TiRT GUI: {url}")
    if not args.no_browser:
        threading.Timer(0.8, webbrowser.open, args=(url,)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping TiRT GUI")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
