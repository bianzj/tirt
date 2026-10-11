"""API smoke tests for source and packaged TIRTEB GUI servers."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def reject_nonfinite_json(value: str):
    raise ValueError(f"non-finite JSON value: {value}")


def request(method: str, url: str, payload: dict | None = None, timeout: int = 45):
    data = None if payload is None else json.dumps(payload).encode("utf-8")
    headers = {} if payload is None else {"Content-Type": "application/json"}
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            body = response.read().decode("utf-8", errors="replace")
            return response.status, json.loads(body, parse_constant=reject_nonfinite_json) if "json" in response.headers.get("Content-Type", "") else body
    except urllib.error.HTTPError as error:
        body = error.read().decode("utf-8", errors="replace")
        try:
            return error.code, json.loads(body, parse_constant=reject_nonfinite_json)
        except json.JSONDecodeError:
            return error.code, body


def wait_ready(base_url: str) -> None:
    deadline = time.time() + 60
    while time.time() < deadline:
        try:
            status, payload = request("GET", f"{base_url}/api/health", timeout=3)
            if status == 200 and payload.get("ok") and "energy_balance" in payload.get("features", []):
                return
        except Exception:
            pass
        time.sleep(0.5)
    raise RuntimeError("TIRTEB GUI did not become ready")


def payload(**overrides: object) -> dict:
    values = {
        "start_time": "214.0139", "end_time": "214.0764", "max_iter": "20",
        "surface_model": "plane", "vegetation_model": "hom",
        "geometry_mode": "0", "vza": "0;20", "vaa": "0;90",
        "spectral_source": "direct", "wavelengths": "0.55;0.86;10.5",
        "observation_time_mode": "all", "observation_quantity": "temperature",
    }
    values.update(overrides)
    return {"overrides": values}


def assert_run(base_url: str, name: str, request_payload: dict) -> dict:
    status, data = request("POST", f"{base_url}/api/run_time", request_payload)
    if status != 200 or not data.get("ok"):
        raise AssertionError(f"{name}: HTTP {status} {data}")
    if not data.get("series") or not data.get("observations"):
        raise AssertionError(f"{name}: no series or observations")
    print(f"PASS {name}: series={data['count']} observations={data['observation_count']}")
    return data


def assert_angle_run(base_url: str, name: str, **overrides: object) -> dict:
    values = {
        "surface_model": "plane", "vegetation_model": "hom",
        "geometry_mode": "3", "include_principal": "1",
        "include_hemisphere": "1", "principal_vza_step": "30",
        "principal_vza_max": "60", "wavelengths": "8;10.5;12",
        "emissivity_leaf": "0.975;0.985;0.98",
        "emissivity_soil": "0.92;0.955;0.965",
    }
    values.update(overrides)
    status, data = request("POST", f"{base_url}/api/run_angle", {"overrides": values})
    if status != 200 or not data.get("ok") or len(data.get("rows", [])) <= 36:
        raise AssertionError(f"{name}: HTTP {status} {data}")
    print(f"PASS {name}: angle_rows={data['count']}")
    return data

def assert_urban_temperatures(data: dict) -> None:
    fields = (
        "Troof_sunlit_C", "Troof_shaded_C", "Twall_sunlit_C", "Twall_shaded_C",
        "Tstreet_sunlit_C", "Tstreet_shaded_C", "Tleaf_sunlit_C", "Tleaf_shaded_C",
    )
    for row in data["series"]:
        for field in fields:
            value = float(row[field])
            if not math.isfinite(value):
                raise AssertionError(f"urban component temperature is not finite: {field}={row[field]}")
    print(f"PASS urban+turbid component temperatures: {len(fields)} fields")


def assert_urban_diagnostics(data: dict) -> None:
    rows = data.get("intermediate", [])
    if not rows:
        raise AssertionError("urban diagnostics are missing")
    row = rows[-1]
    values = {key: float(row[key]) for key in ("raa_roof", "raa_wall", "raa_street", "raa_leaf")}
    if not (values["raa_roof"] < values["raa_street"] == values["raa_wall"] < values["raa_leaf"]):
        raise AssertionError(f"unexpected urban resistance ordering: {values}")
    print(f"PASS urban GUI diagnostics: {values}")


def assert_urban_daylight(data: dict) -> None:
    wall_delta = max(abs(float(row["Twall_sunlit_C"]) - float(row["Twall_shaded_C"])) for row in data["series"])
    street_delta = max(abs(float(row["Tstreet_sunlit_C"]) - float(row["Tstreet_shaded_C"])) for row in data["series"])
    if wall_delta <= 0.05 or street_delta <= 0.05:
        raise AssertionError(f"urban sunlit/shaded separation missing: wall={wall_delta}, street={street_delta}")
    print(f"PASS urban daylight separation: wall_delta={wall_delta:.3f}; street_delta={street_delta:.3f}")

def run_tests(base_url: str) -> None:
    wait_ready(base_url)
    status, html = request("GET", f"{base_url}/time.html")
    assert status == 200 and "Energy Balance" in html
    print("PASS index page")
    assert_angle_run(base_url, "angle plane+hom")
    assert_angle_run(base_url, "angle urban+bare", surface_model="urban", vegetation_model="bare")
    assert_run(base_url, "default", payload())
    for surface, vegetation in (("plane", "bare"), ("plane", "hom"), ("plane", "row"), ("plane", "crown"), ("slope", "hom"), ("terrain", "bare"), ("terrain", "crown"), ("urban", "bare"), ("urban", "hom"), ("urban", "crown")):
        data = assert_run(base_url, f"{surface}+{vegetation}", payload(surface_model=surface, vegetation_model=vegetation))
        if surface == "urban" and vegetation == "hom":
            assert_urban_temperatures(data)
            assert_urban_diagnostics(data)
            daylight = assert_run(base_url, "urban+hom daylight", payload(surface_model=surface, vegetation_model=vegetation, start_time="214.4514", end_time="214.4514"))
            assert_urban_daylight(daylight)

    files = {
        "directions": "0 0\n20 90\n40 180\n",
        "wavelengths": "0.55\n0.86\n10.5\n",
        "leaf_optics": "0.45 0.10 0.05\n0.86 0.45 0.40\n2.40 0.05 0.03\n",
        "soil_optics": "0.45 0.20\n0.86 0.30\n2.40 0.35\n",
        "urban_spectrum": "0.45 0.05 0.08 0.045 0.95 0.92 0.955\n10.5 0.05 0.08 0.045 0.95 0.92 0.955\n12.0 0.05 0.08 0.045 0.95 0.92 0.955\n",
    }
    status, data = request("POST", f"{base_url}/api/run_time", {**payload(surface_model="urban", vegetation_model="hom"), "files": files})
    if status != 200 or not data.get("ok") or data.get("observation_count", 0) != data.get("count", 0) * 3:
        raise AssertionError(f"uploaded files: HTTP {status} {data}")
    print(f"PASS uploaded files: series={data['count']} observations={data['observation_count']}")

    status, data = request("POST", f"{base_url}/api/run_time", payload(observation_time_mode="time", observation_times="214.0139"))
    if status != 200 or data.get("observation_count") != 4:
        raise AssertionError(f"selected time: HTTP {status} {data}")
    print("PASS selected observation time")


def start_server(mode: str, port: int) -> subprocess.Popen:
    if mode == "exe":
        command = [str(ROOT / "dist" / "TiRT-service.exe"), "--server-only", "--port", str(port)]
    else:
        command = [sys.executable, str(ROOT / "gui" / "server.py"), "--host", "127.0.0.1", "--port", str(port)]
    return subprocess.Popen(command, cwd=ROOT, stdout=subprocess.DEVNULL)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("source", "exe"), default="source")
    parser.add_argument("--port", type=int, default=8899)
    args = parser.parse_args()
    process = start_server(args.mode, args.port)
    try:
        run_tests(f"http://127.0.0.1:{args.port}")
    finally:
        process.terminate()
        try:
            process.wait(timeout=8)
        except subprocess.TimeoutExpired:
            process.kill()


if __name__ == "__main__":
    main()
