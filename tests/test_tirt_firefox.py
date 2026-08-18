from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait

ROOT = Path(__file__).resolve().parent.parent
PYTHON = Path(r"C:\work\miniconda\envs\python311\python.exe")
FIREFOX = Path(r"C:\Program Files\Mozilla Firefox\firefox.exe")


def start_server(port: int) -> subprocess.Popen:
    return subprocess.Popen(
        [str(PYTHON), str(ROOT / "gui" / "server.py"), "--host", "127.0.0.1", "--port", str(port)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Test TiRT GUI with installed Firefox")
    parser.add_argument("--port", type=int, default=8785)
    parser.add_argument("--screenshot", type=Path, default=ROOT / "docs" / "firefox_results.png")
    args = parser.parse_args()

    proc = start_server(args.port)
    driver = None
    try:
        time.sleep(1.0)
        options = Options()
        options.binary_location = str(FIREFOX)
        options.add_argument("-headless")
        driver = webdriver.Firefox(options=options)
        driver.set_window_size(1440, 1000)
        driver.get(f"http://127.0.0.1:{args.port}/")
        wait = WebDriverWait(driver, 60)
        wait.until(EC.presence_of_element_located((By.ID, "run-button")))
        assert "TiRT Directional Simulator" in driver.title
        driver.find_element(By.ID, "run-button").click()
        wait.until(lambda d: "is-hidden" not in d.find_element(By.ID, "results-content").get_attribute("class"))
        wait.until(lambda d: d.find_element(By.ID, "observation-count").text != "--")
        status = driver.find_element(By.CSS_SELECTOR, "#status-pill span").text
        count = driver.find_element(By.ID, "observation-count").text
        mean = driver.find_element(By.ID, "mean-bt").text
        assert count and count != "--"
        args.screenshot.parent.mkdir(parents=True, exist_ok=True)
        driver.save_screenshot(str(args.screenshot))
        print(f"PASS installed Firefox status={status} observations={count} mean={mean} screenshot={args.screenshot}")
    finally:
        if driver is not None:
            driver.quit()
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()


if __name__ == "__main__":
    main()
