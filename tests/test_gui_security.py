"""The local GUI API must reject foreign paths, hosts and origins."""
import http.client
import json
import threading
from http.server import ThreadingHTTPServer

import pytest

from gui.server import Handler, PROJECT_ROOT
from gui.tirteb_server import _check_paths


@pytest.fixture(scope="module")
def server():
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield httpd.server_address[1]
    httpd.shutdown()
    httpd.server_close()


def request(port, method, path, headers=None, body=None):
    headers = dict(headers or {})
    headers.setdefault("Host", f"127.0.0.1:{port}")
    data = None if body is None else json.dumps(body).encode("utf-8")
    if data is not None:
        headers["Content-Type"] = "application/json"
    connection = http.client.HTTPConnection("127.0.0.1", port, timeout=60)
    try:
        connection.request(method, path, body=data, headers=headers)
        response = connection.getresponse()
        return response.status, response.read()
    finally:
        connection.close()


def test_data_and_work_dir_paths_are_accepted(tmp_path):
    _check_paths({
        "spectrals_file": "data/spectrals/optical_properties.txt",
        "hemisphere_geometry_file": str(PROJECT_ROOT / "data" / "directions" / "hemisphere_directions.txt"),
        "geometry_file": str(tmp_path / "geometry.txt"),
        "meteo": str(tmp_path / "meteo.txt"),
        "additional_geometry_file": "",
        "lai": "2",
    }, tmp_path)


@pytest.mark.parametrize("key, value", [
    ("spectrals_file", r"C:\Windows\win.ini"),
    ("spectrals_file", "data/../run.py"),
    ("geometry_file", "../outside.txt"),
    ("meteo", r"C:\Windows\win.ini"),
    ("result_dir", r"C:\Temp"),
])
def test_paths_outside_data_are_rejected(tmp_path, key, value):
    with pytest.raises(ValueError, match=key):
        _check_paths({key: value}, tmp_path)


@pytest.mark.parametrize("path, overrides", [
    ("/api/export", {"spectral_source": "file", "spectrals_file": r"C:\Windows\win.ini"}),
    ("/api/run", {"geometry_mode": 2, "geometry_file": r"C:\Windows\win.ini"}),
    ("/api/run_time", {"meteo": r"C:\Windows\win.ini"}),
])
def test_api_rejects_foreign_paths(server, path, overrides):
    status, body = request(server, "POST", path, body={"overrides": overrides})
    assert status == 400
    payload = json.loads(body)
    assert payload["ok"] is False and "data folder" in payload["error"]
    assert b"[fonts]" not in body


@pytest.mark.parametrize("host", ["127.0.0.1", "localhost", "LOCALHOST"])
def test_health_with_local_host(server, host):
    status, body = request(server, "GET", "/api/health", {"Host": f"{host}:{server}"})
    assert status == 200 and json.loads(body)["ok"] is True


@pytest.mark.parametrize("host", ["evil.example:{port}", "127.0.0.1:1", "127.0.0.1", "attacker.localhost:{port}"])
def test_foreign_host_is_forbidden(server, host):
    status, _ = request(server, "GET", "/api/health", {"Host": host.format(port=server)})
    assert status == 403


def test_same_origin_is_allowed(server):
    status, _ = request(server, "POST", "/api/run", {"Origin": f"http://127.0.0.1:{server}"},
                        body={"overrides": {"spectrals_file": r"C:\Windows\win.ini"}})
    assert status == 400


@pytest.mark.parametrize("origin", ["http://evil.example", "null", "https://127.0.0.1:{port}", "http://localhost:1"])
def test_cross_origin_is_forbidden(server, origin):
    origin = origin.format(port=server)
    assert request(server, "GET", "/api/health", {"Origin": origin})[0] == 403
    assert request(server, "POST", "/api/export", {"Origin": origin}, body={"overrides": {}})[0] == 403


def test_cross_origin_shutdown_is_rejected(server):
    status, _ = request(server, "POST", "/api/shutdown", {"Origin": "http://evil.example"})
    assert status == 403
    status, _ = request(server, "POST", "/api/shutdown", {"Host": "evil.example"})
    assert status == 403
    assert request(server, "GET", "/api/health")[0] == 200


def test_time_run_input_keeps_crown_names_and_drops_tirt_geometry(tmp_path, monkeypatch):
    import csv
    from gui import tirteb_server

    captured, appended = {}, set()

    def fake_run(input_path, output_path):
        with open(input_path, newline="", encoding="utf-8") as file:
            for row in csv.reader(file):
                if len(row) >= 3:
                    captured[row[1]] = row[2]
                    if row[0] == "gui":
                        appended.add(row[1])
        raise RuntimeError("stop")

    monkeypatch.setattr(tirteb_server, "run", fake_run)
    overrides = {"crown_rad_a": "2.5", "crown_rad_b": "3", "crown_height": "4", "slope_aspect": "90",
                 "geometry_mode": "0", "geometry_pair_mode": "zip", "vza": "0;30", "vaa": "0;90",
                 "principal_vza_step": "10", "additional_geometry_source": "none"}
    with pytest.raises(RuntimeError, match="stop"):
        tirteb_server._time_request({"overrides": overrides}, tmp_path)
    for key in ("geometry_mode", "geometry_pair_mode", "principal_vza_step", "additional_geometry_source",
                "crown_radius", "crown_width", "slope_aspect"):
        assert key not in appended
    assert (captured["crown_rad_a"], captured["crown_rad_b"], captured["crown_height"]) == ("2.5", "3", "4")
    assert captured["terrain_aspect"] == "90"
