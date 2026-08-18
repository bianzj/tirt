from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parent.parent
PYTHON = Path(r"C:\work\miniconda\envs\python311\python.exe")


def start_server(port: int) -> subprocess.Popen:
    return subprocess.Popen(
        [str(PYTHON), str(ROOT / "gui" / "server.py"), "--host", "127.0.0.1", "--port", str(port)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Browser-level TiRT GUI test")
    parser.add_argument("--browser", choices=("chromium", "firefox"), default="firefox")
    parser.add_argument("--port", type=int, default=8784)
    parser.add_argument("--screenshot", type=Path, default=ROOT / "docs" / "browser_firefox_results.png")
    args = parser.parse_args()

    proc = start_server(args.port)
    errors: list[str] = []
    try:
        with sync_playwright() as p:
            browser_type = getattr(p, args.browser)
            browser = browser_type.launch(headless=True)
            page = browser.new_page(viewport={"width": 1440, "height": 1000})
            page.on("console", lambda msg: errors.append(f"console {msg.type}: {msg.text}") if msg.type == "error" else None)
            page.on("pageerror", lambda exc: errors.append(f"pageerror: {exc}"))
            page.goto(f"http://127.0.0.1:{args.port}/", wait_until="networkidle", timeout=45000)
            page.wait_for_selector("#run-button", timeout=30000)
            page.click("#run-button")
            page.wait_for_selector("#results-content:not(.is-hidden)", timeout=60000)
            page.wait_for_function("() => document.querySelector('#observation-count')?.textContent !== '--'", timeout=30000)
            count = page.locator("#observation-count").inner_text()
            mean = page.locator("#mean-bt").inner_text()
            status = page.locator("#status-pill span").inner_text()
            page.screenshot(path=str(args.screenshot), full_page=True)
            browser.close()
        if errors:
            raise RuntimeError("; ".join(errors))
        print(f"PASS browser={args.browser} status={status} observations={count} mean={mean} screenshot={args.screenshot}")
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()


if __name__ == "__main__":
    main()
