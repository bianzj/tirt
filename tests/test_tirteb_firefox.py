"""Run the temporal GUI through the installed Firefox browser."""

from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait

ROOT = Path(__file__).resolve().parents[1]
PYTHON = Path(r"C:\work\miniconda\envs\python311\python.exe")
FIREFOX = Path(r"C:\Program Files\Mozilla Firefox\firefox.exe")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8902)
    parser.add_argument("--screenshot", type=Path, default=ROOT / "docs" / "screenshots" / "firefox_results.png")
    parser.add_argument("--surface", choices=("plane", "urban"), default="plane")
    args = parser.parse_args()
    server = subprocess.Popen([str(PYTHON), str(ROOT / "gui" / "server.py"), "--host", "127.0.0.1", "--port", str(args.port)], cwd=ROOT)
    driver = None
    try:
        time.sleep(1.0)
        options = Options()
        options.binary_location = str(FIREFOX)
        options.add_argument("-headless")
        driver = webdriver.Firefox(options=options)
        driver.set_window_size(1440, 1000)
        driver.get(f"http://127.0.0.1:{args.port}/time.html")
        wait = WebDriverWait(driver, 90)
        wait.until(lambda item: item.find_element(By.ID, "run-button").is_displayed())
        assert "STREAM TIRT Temporal Simulator" in driver.title
        summaries = driver.find_elements(By.CSS_SELECTOR, "details.form-section > summary")
        summaries[2].click()
        Select(driver.find_element(By.ID, "surface-model")).select_by_value(args.surface)
        Select(driver.find_element(By.ID, "vegetation-model")).select_by_value("hom")
        summaries[0].click()
        time_range = driver.find_element(By.ID, "time-range")
        time_range.clear()
        time_range.send_keys("214.4514-214.4514")
        driver.find_element(By.ID, "run-button").click()
        wait.until(lambda item: item.find_element(By.ID, "status-pill").text in {"Complete", "Error"})
        status_element = driver.find_element(By.ID, "status-pill")
        status = status_element.text
        error_title = status_element.get_attribute("title") or ""
        toast_text = driver.find_element(By.ID, "toast").get_attribute("textContent") or ""
        args.screenshot.parent.mkdir(parents=True, exist_ok=True)
        driver.save_screenshot(str(args.screenshot))
        print(f"Firefox status={status!r} title={error_title!r} toast={toast_text!r}")
        assert status == "Complete", error_title or toast_text
        assert driver.find_element(By.ID, "diagnostics-card").is_displayed() == (args.surface == "urban")
        count = driver.find_element(By.ID, "metric-count").text
        assert count != "--" and count != "0"
        print(f"PASS Firefox {args.surface}+turbid: status=Complete time_steps={count} screenshot={args.screenshot}")
    finally:
        if driver is not None:
            driver.quit()
        server.terminate()
        try:
            server.wait(timeout=8)
        except subprocess.TimeoutExpired:
            server.kill()


if __name__ == "__main__":
    main()
