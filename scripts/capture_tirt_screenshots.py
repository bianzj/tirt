from __future__ import annotations

import subprocess
import time
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select, WebDriverWait

ROOT = Path(__file__).resolve().parent.parent
PYTHON = Path(r"C:\work\miniconda\envs\python311\python.exe")
FIREFOX = Path(r"C:\Program Files\Mozilla Firefox\firefox.exe")
OUT = ROOT / "docs" / "screenshots"


def start_server(port: int) -> subprocess.Popen:
    return subprocess.Popen(
        [str(PYTHON), str(ROOT / "gui" / "server.py"), "--host", "127.0.0.1", "--port", str(port)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def screenshot(driver, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    driver.save_screenshot(str(OUT / name))


def run_current(driver, wait: WebDriverWait) -> None:
    driver.find_element(By.ID, "run-button").click()
    wait.until(lambda d: "is-hidden" not in d.find_element(By.ID, "results-content").get_attribute("class"))
    wait.until(lambda d: d.find_element(By.ID, "observation-count").text != "--")
    time.sleep(0.4)


def main() -> None:
    port = 8786
    proc = start_server(port)
    driver = None
    try:
        time.sleep(1.0)
        options = Options()
        options.binary_location = str(FIREFOX)
        options.add_argument("-headless")
        driver = webdriver.Firefox(options=options)
        wait = WebDriverWait(driver, 60)

        driver.set_window_size(1440, 1000)
        driver.get(f"http://127.0.0.1:{port}/")
        wait.until(EC.presence_of_element_located((By.ID, "run-button")))
        screenshot(driver, "01_home.png")

        run_current(driver, wait)
        screenshot(driver, "02_default_results.png")

        Select(driver.find_element(By.ID, "surface-model")).select_by_value("urban")
        Select(driver.find_element(By.ID, "vegetation-model")).select_by_value("crown")
        driver.execute_script("document.querySelector('#building-height').value='35'; document.querySelector('#crown-density').value='0.04';")
        screenshot(driver, "03_urban_crown_settings.png")
        run_current(driver, wait)
        screenshot(driver, "04_urban_crown_results.png")

        driver.execute_script("document.querySelector('input[name=additional-source][value=manual]').click(); document.querySelector('#include-hemisphere').checked=false; document.querySelector('#include-principal').checked=false; document.querySelector('#manual-geometry').value='0/0; 20/90; 40/180; 60/270';")
        run_current(driver, wait)
        screenshot(driver, "05_manual_geometry_results.png")

        driver.set_window_size(390, 900)
        time.sleep(0.5)
        screenshot(driver, "06_mobile_view.png")
        print(f"screenshots={OUT}")
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
