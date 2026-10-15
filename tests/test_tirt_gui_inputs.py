"""Isolated API and portable-input regression tests; all outputs stay temporary."""
from __future__ import annotations

import csv
import io
import json
import shutil
import subprocess
import sys
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
import zipfile
from http.server import ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from base.inputio import read_input, resolve_path  # noqa: E402
from gui.server import Handler  # noqa: E402
from run import run  # noqa: E402


class GuiInputsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.url = f"http://127.0.0.1:{cls.server.server_port}"

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join(timeout=5)

    def request(self, endpoint, payload):
        request = urllib.request.Request(self.url + endpoint, data=json.dumps(payload).encode("utf-8"),
                                         headers={"Content-Type": "application/json"})
        try:
            response = urllib.request.urlopen(request, timeout=45)
        except urllib.error.HTTPError as error:
            response = error
        with response:
            body = response.read()
            return response.status, response.headers, (
                json.loads(body) if response.headers.get_content_type() == "application/json" else body
            )

    @staticmethod
    def overrides(**extra):
        return {"surface_model": "plane", "vegetation_model": "bare", "geometry_mode": "3",
                "include_principal": "0", "include_hemisphere": "0", "principal_vza_max": "10",
                "additional_geometry_source": "manual", "observation_vza": "15", "observation_vaa": "77", **extra}

    def test_manual_selection_ignores_retained_uploaded_geometry(self):
        status, _, result = self.request("/api/run", {"overrides": self.overrides(), "geometry_text": "40 180\n"})
        self.assertEqual(status, 200, result)
        self.assertEqual({(row["vza"], row["vaa"]) for row in result["rows"]}, {(15, 77)})

    def test_hemisphere_and_uploaded_directions_are_both_kept(self):
        status, _, result = self.request("/api/run", {"overrides": self.overrides(
            include_hemisphere="1", additional_geometry_source="file"), "geometry_text": "45 77\n"})
        self.assertEqual(status, 200, result)
        directions = {(row["vza"], row["vaa"]) for row in result["rows"]}
        self.assertIn((10, 30), directions)
        self.assertIn((45, 77), directions)
        self.assertEqual(len(directions), 14)

    def test_hemisphere_and_manual_ignore_retained_upload(self):
        status, _, result = self.request("/api/run", {"overrides": self.overrides(include_hemisphere="1"),
                                                    "geometry_text": "40 180\n"})
        self.assertEqual(status, 200, result)
        directions = {(row["vza"], row["vaa"]) for row in result["rows"]}
        self.assertIn((10, 30), directions)
        self.assertIn((15, 77), directions)
        self.assertNotIn((40, 180), directions)

    def test_export_replays_after_all_upload_temporary_files_are_deleted(self):
        payload = {"overrides": self.overrides(additional_geometry_source="file", include_hemisphere="1"),
                   "geometry_text": "45 77 20 80\n", "files": {
                       "wavelengths": "8\n10.5\n12\n",
                       "leaf_spectrum": "8 0.83\n10.5 0.87\n12 0.89\n",
                       "soil_spectrum": "8 0.71\n10.5 0.77\n12 0.79\n",
                       "terrain_spectrum": "8 0.75\n10.5 0.8\n12 0.9\n",
                       "roof_spectrum": "8 0.8\n10.5 0.85\n12 0.9\n",
                       "wall_spectrum": "8 0.79\n10.5 0.86\n12 0.91\n",
                       "street_spectrum": "8 0.82\n10.5 0.84\n12 0.87\n",
                   }}
        status, _, result = self.request("/api/run", payload)
        self.assertEqual(status, 200, result)
        status, headers, bundle = self.request("/api/export", payload)
        self.assertEqual(status, 200, bundle)
        self.assertEqual(headers.get_content_type(), "application/zip")
        with tempfile.TemporaryDirectory(prefix="tirt-portable-input-") as directory:
            extracted = Path(directory) / "relocated" / "case"
            with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
                self.assertIn("input.csv", archive.namelist())
                self.assertIn("README.txt", archive.namelist())
                archive.extractall(extracted)
            config = read_input(extracted / "input.csv")
            for key in ("spectrals_file", "vegetation_spectrum_file", "soil_spectrum_file", "terrain_spectrum_file",
                        "roof_spectrum_file", "wall_spectrum_file", "street_spectrum_file",
                        "hemisphere_geometry_file", "additional_geometry_file"):
                self.assertFalse(Path(config[key]).is_absolute(), key)
                self.assertTrue(resolve_path(extracted / "input.csv", config[key]).is_file(), key)
            output = run(extracted / "input.csv", extracted / "replay" / "output.csv")
            with output.open(newline="", encoding="utf-8") as stream:
                replay = list(csv.DictReader(stream))
            self.assertEqual(len(replay), len(result["rows"]))
            for actual, expected in zip(replay, result["rows"]):
                for key in ("vza", "vaa", "sza", "saa", "wavelength_um", "radiance", "brightness_temperature_K"):
                    self.assertAlmostEqual(float(actual[key]), float(expected[key]), places=12, msg=key)

    def test_malformed_request_is_reported(self):
        for payload in ([], {"overrides": []}, {"files": []}, {"files": {"soil_spectrum": [1, 2]}}):
            status, _, result = self.request("/api/run", payload)
            self.assertEqual(status, 400, result)
            self.assertFalse(result["ok"])

    def test_frontend_input_regressions(self):
        node = shutil.which("node")
        if node is None:
            self.skipTest("Node.js is required for the GUI JavaScript checks")
        result = subprocess.run([node, str(ROOT / "tests" / "test_tirt_gui_inputs.cjs")], cwd=ROOT,
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
