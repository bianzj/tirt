"""Capture reproducible TIRTEB GUI states with installed Firefox."""

from __future__ import annotations

import subprocess
import time
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.ui import WebDriverWait

ROOT = Path(__file__).resolve().parents[1]
PYTHON = Path(r"C:\work\miniconda\envs\python311\python.exe")
FIREFOX = Path(r"C:\Program Files\Mozilla Firefox\firefox.exe")
PORT = 8905
OUT = ROOT / "docs" / "screenshots"


def save(driver: webdriver.Firefox, name: str) -> None:
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    driver.save_screenshot(str(path))
    print(path)


def run_and_wait(driver: webdriver.Firefox, wait: WebDriverWait) -> None:
    driver.find_element(By.ID, "run-button").click()
    wait.until(lambda item: item.find_element(By.ID, "status-pill").text == "Complete")
    time.sleep(0.5)


def select_value(driver: webdriver.Firefox, element_id: str, value: str) -> None:
    """Select values without requiring the option to be visible in headless mode."""
    driver.execute_script(
        """
        const select = document.getElementById(arguments[0]);
        select.value = arguments[1];
        select.dispatchEvent(new Event('change', {bubbles: true}));
        """,
        element_id,
        value,
    )


def set_value(driver: webdriver.Firefox, element_id: str, value: str) -> None:
    """Set a form value without requiring the control to be visible."""
    driver.execute_script(
        """
        const input = document.getElementById(arguments[0]);
        input.value = arguments[1];
        input.dispatchEvent(new Event('input', {bubbles: true}));
        input.dispatchEvent(new Event('change', {bubbles: true}));
        """,
        element_id,
        value,
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    server = subprocess.Popen([str(PYTHON), str(ROOT / "gui" / "server.py"), "--host", "127.0.0.1", "--port", str(PORT)], cwd=ROOT)
    driver = None
    try:
        time.sleep(1.0)
        options = Options()
        options.binary_location = str(FIREFOX)
        options.add_argument("-headless")
        driver = webdriver.Firefox(options=options)
        driver.set_window_size(1440, 1000)
        driver.get(f"http://127.0.0.1:{PORT}/time.html")
        wait = WebDriverWait(driver, 120)
        wait.until(lambda item: item.find_element(By.ID, "run-button").is_displayed())
        save(driver, "07_energy_balance.png")

        set_value(driver, "time-range", "214.0139-214.0764")
        run_and_wait(driver, wait)
        save(driver, "08_energy_balance_results.png")

        select_value(driver, "surface-model", "urban")
        select_value(driver, "vegetation-model", "crown")
        time.sleep(0.5)
        save(driver, "09_energy_balance_urban_settings.png")
        run_and_wait(driver, wait)
        save(driver, "10_energy_balance_urban_results.png")

        Select(driver.find_element(By.ID, "surface-model")).select_by_value("terrain")
        set_value(driver, "observation-angles", "0/0; 20/90; 40/180")
        set_value(driver, "wavelengths", "0.55; 0.86; 10.5")
        run_and_wait(driver, wait)
        save(driver, "11_energy_balance_manual_geometry.png")

        driver.set_window_size(390, 844)
        save(driver, "12_energy_balance_mobile.png")
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
