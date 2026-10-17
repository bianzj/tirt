"""Regress the launcher server, actual temporal solver, and GUI parameter contract."""
import json
import math
import shutil
import subprocess
import sys
import threading
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gui.server import Handler
from gui.tirteb_server import _time_overrides
from scripts.check_dynamic_api import check_dynamic_api


class DynamicGuiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.url = f'http://127.0.0.1:{cls.server.server_port}'
        cls.opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join(timeout=5)

    def request(self, path, body=None):
        req = urllib.request.Request(self.url+path, data=None if body is None else json.dumps(body).encode(),
                                     headers={'Content-Type': 'application/json'})
        try:
            response = self.opener.open(req, timeout=180)
        except urllib.error.HTTPError as error:
            if path == '/api/run_time' and isinstance(body, dict) and body.get('meteo_text'):
                raise AssertionError(error.read().decode()) from error
            raise
        with response:
            return response.headers, response.read()

    def test_time_calculations_and_uploaded_files(self):
        result = check_dynamic_api(self.request)
        self.assertEqual(len(result['cases']), 6)

    def test_invalid_time_requests_return_json_400(self):
        for payload in ([], {'overrides': []}, {'files': ['bad']},
                        {'overrides': {'emissivity_leaf': 1.2}},
                        {'overrides': {'geometry_pair_mode': 'zip', 'vza': '20;30', 'vaa': '15'}}):
            with self.assertRaises(urllib.error.HTTPError) as result:
                self.request('/api/run_time', payload)
            self.assertEqual(result.exception.code, 400)
            self.assertFalse(json.loads(result.exception.read())['ok'])

    def test_gui_values_reach_solver_names(self):
        result = _time_overrides({'overrides': dict(canopy_height=7, crown_rad_a=2, crown_rad_b=4,
                crown_height=3, slope_angle=25, slope_aspect=120, emissivity_leaf=.98,
                emissivity_soil=.95, emissivity_roof=.9)})
        for key, expected in dict(height=7, crown_width=2, crown_radius=2, crown_height=4,
                                  terrain_slope=25, terrain_aspect=120, urban_emissivity_roof=.9).items():
            self.assertEqual(result[key], expected)
        self.assertAlmostEqual(float(result['leaf_reflectance']), .02)
        self.assertAlmostEqual(float(result['soil_reflectance']), .05)

    def test_legacy_angle_endpoint_remains_available(self):
        result = json.loads(self.request('/api/run_angle', {'overrides': dict(surface_model='plane',
            vegetation_model='bare', geometry_mode=0, vza=0, vaa=0, wavelengths='10.5')})[1])
        self.assertTrue(result['ok'])
        self.assertGreater(result['count'], 0)

    def test_frontend_identity_and_payload(self):
        node = shutil.which('node')
        if not node:
            self.skipTest('Node.js is needed for frontend checks')
        result = subprocess.run([node, str(ROOT/'tests/test_tirt_dynamic_gui.cjs')],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
