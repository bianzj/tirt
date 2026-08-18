from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PYTHON = Path(r"C:\work\miniconda\envs\python311\python.exe")


@dataclass
class Result:
    name: str
    ok: bool
    detail: str


def request(method: str, url: str, payload: dict | None = None, timeout: int = 30):
    data = None if payload is None else json.dumps(payload).encode("utf-8")
    headers = {} if payload is None else {"Content-Type": "application/json"}
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read()
            text = body.decode("utf-8", errors="replace")
            if resp.headers.get("Content-Type", "").startswith("application/json"):
                return resp.status, json.loads(text)
            return resp.status, text
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        try:
            return exc.code, json.loads(body)
        except json.JSONDecodeError:
            return exc.code, body


def wait_ready(base_url: str, seconds: int = 45) -> dict:
    last = None
    deadline = time.time() + seconds
    while time.time() < deadline:
        try:
            status, data = request("GET", f"{base_url}/api/health", timeout=3)
            if status == 200 and isinstance(data, dict) and data.get("ok"):
                return data
            last = f"status={status} data={data}"
        except Exception as exc:  # noqa: BLE001
            last = str(exc)
        time.sleep(0.5)
    raise RuntimeError(f"server not ready: {last}")


def add(results: list[Result], name: str, condition: bool, detail: str) -> None:
    results.append(Result(name, condition, detail))


def run_api_tests(base_url: str) -> list[Result]:
    results: list[Result] = []
    health = wait_ready(base_url)
    add(results, "health endpoint", health.get("ok") is True and health.get("engine") == "tirt", str(health))

    status, html = request("GET", f"{base_url}/")
    add(results, "index page", status == 200 and "TiRT Directional Simulator" in str(html), f"HTTP {status}")

    for asset in ("app.js", "styles.css", "vendor/three.module.js", "vendor/controls/OrbitControls.js"):
        status, body = request("GET", f"{base_url}/{asset}")
        add(results, f"static asset {asset}", status == 200 and len(str(body)) > 100, f"HTTP {status}, bytes~{len(str(body))}")

    required_fields = {"vza", "vaa", "sza", "saa", "raa", "wavelength_um", "radiance", "brightness_temperature_K", "brightness_temperature_C"}

    status, data = request("POST", f"{base_url}/api/run", {})
    rows = data.get("rows", []) if isinstance(data, dict) else []
    add(results, "default run", status == 200 and data.get("ok") and data.get("count") == len(rows) and len(rows) > 0, f"HTTP {status}, count={data.get('count') if isinstance(data, dict) else None}")
    add(results, "output fields", bool(rows) and required_fields.issubset(rows[0]), f"fields={sorted(rows[0]) if rows else []}")

    matrix_cases = [
        ("plane", "bare"), ("plane", "hom"), ("plane", "row"), ("plane", "crown"),
        ("slope", "hom"), ("terrain", "bare"), ("terrain", "crown"), ("urban", "bare"), ("urban", "crown"),
    ]
    for surface, vegetation in matrix_cases:
        overrides = {
            "surface_model": surface,
            "vegetation_model": vegetation,
            "geometry_mode": "0",
            "vza": "0;30",
            "vaa": "0;180",
        }
        status, data = request("POST", f"{base_url}/api/run", {"overrides": overrides})
        add(results, f"model {surface}+{vegetation}", status == 200 and data.get("ok") and data.get("count", 0) > 0, f"HTTP {status}, count={data.get('count') if isinstance(data, dict) else None}")

    manual = {
        "overrides": {
            "geometry_mode": "3",
            "include_principal": "0",
            "include_hemisphere": "0",
            "additional_geometry_source": "manual",
            "observation_vza": "0;20;40",
            "observation_vaa": "0;90;180",
            "wavelengths": "10.5",
            "emissivity_leaf": "0.985",
            "emissivity_soil": "0.955",
            "emissivity_terrain": "0.950",
            "emissivity_roof": "0.950",
            "emissivity_wall": "0.930",
            "emissivity_street": "0.955",
        }
    }
    status, data = request("POST", f"{base_url}/api/run", manual)
    add(results, "manual geometry", status == 200 and data.get("ok") and data.get("count") == 3, f"HTTP {status}, count={data.get('count') if isinstance(data, dict) else None}")

    upload_geometry = {
        "overrides": {"geometry_mode": "2", "wavelengths": "10.5", "emissivity_leaf": "0.985", "emissivity_soil": "0.955", "emissivity_terrain": "0.950", "emissivity_roof": "0.950", "emissivity_wall": "0.930", "emissivity_street": "0.955"},
        "geometry_text": "0 0 30 0\n25 90 30 0\n50 180 30 0\n",
    }
    status, data = request("POST", f"{base_url}/api/run", upload_geometry)
    add(results, "uploaded geometry text", status == 200 and data.get("ok") and data.get("count") == 3, f"HTTP {status}, count={data.get('count') if isinstance(data, dict) else None}")

    spectra_files = {
        "overrides": {"geometry_mode": "0", "vza": "0", "vaa": "0"},
        "files": {
            "wavelengths": "8\n10.5\n12\n",
            "leaf_spectrum": "8 0.975\n10.5 0.985\n12 0.980\n",
            "soil_spectrum": "8 0.920\n10.5 0.955\n12 0.965\n",
            "roof_spectrum": "8 0.920\n10.5 0.950\n12 0.940\n",
            "wall_spectrum": "8 0.900\n10.5 0.930\n12 0.920\n",
            "street_spectrum": "8 0.940\n10.5 0.955\n12 0.950\n",
        },
    }
    status, data = request("POST", f"{base_url}/api/run", spectra_files)
    add(results, "uploaded spectral files", status == 200 and data.get("ok") and data.get("count", 0) >= 3, f"HTTP {status}, count={data.get('count') if isinstance(data, dict) else None}")

    status, data = request("POST", f"{base_url}/api/run", {"overrides": {"surface_model": "bad"}})
    add(results, "invalid input error", status == 400 and data.get("ok") is False and "surface_model" in data.get("error", ""), f"HTTP {status}, error={data.get('error') if isinstance(data, dict) else None}")

    return results


def start_server(mode: str, port: int) -> subprocess.Popen:
    if mode == "exe":
        cmd = [str(ROOT / "dist" / "TiRT.exe"), "--no-browser", "--port", str(port)]
    else:
        python = str(PYTHON if PYTHON.exists() else sys.executable)
        cmd = [python, str(ROOT / "gui" / "server.py"), "--host", "127.0.0.1", "--port", str(port)]
    return subprocess.Popen(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)


def write_report(mode: str, port: int, results: list[Result], report_path: Path) -> None:
    passed = sum(1 for item in results if item.ok)
    lines = [
        "# TiRT GUI 测试报告",
        "",
        f"- 测试模式：`{mode}`",
        f"- 测试地址：`http://127.0.0.1:{port}`",
        f"- 通过：`{passed}/{len(results)}`",
        "",
        "| 测试项 | 结果 | 说明 |",
        "|---|---:|---|",
    ]
    for item in results:
        lines.append(f"| {item.name} | {'PASS' if item.ok else 'FAIL'} | `{item.detail.replace('|', '/')}` |")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Test TiRT GUI HTTP/API behavior")
    parser.add_argument("--mode", choices=("source", "exe"), default="source")
    parser.add_argument("--port", type=int, default=8780)
    parser.add_argument("--report", type=Path, default=ROOT / "docs" / "gui_test_report.md")
    args = parser.parse_args()

    proc = start_server(args.mode, args.port)
    try:
        results = run_api_tests(f"http://127.0.0.1:{args.port}")
        write_report(args.mode, args.port, results, args.report)
        for item in results:
            print(f"{'PASS' if item.ok else 'FAIL'} {item.name}: {item.detail}")
        if not all(item.ok for item in results):
            raise SystemExit(1)
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()


if __name__ == "__main__":
    main()
