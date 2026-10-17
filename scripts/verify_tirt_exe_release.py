"""Test the real windowed TiRT executable with no Python/Conda paths inherited."""
import argparse
import csv
import hashlib
import io
import json
import math
import os
import socket
import struct
import subprocess
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    exe, bundle, output = args.exe.resolve(), args.bundle.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    for key in list(env):
        if key.startswith(('PYTHON', 'CONDA', '_PYI')):
            env.pop(key, None)
    windows = Path(env.get('SystemRoot', 'C:/Windows'))
    env['PATH'] = os.pathsep.join(map(str, [windows/'System32', windows, windows/'System32/Wbem']))
    startup = dict(cwd=output, env=env, creationflags=subprocess.CREATE_NO_WINDOW)
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    report = dict(exe=str(exe), sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                  executable_size=exe.stat().st_size, python_and_conda_paths_removed=True)
    with exe.open('rb') as f:
        f.seek(0x3c)
        pe = struct.unpack('<I', f.read(4))[0]
        f.seek(pe+4)
        report['pe_machine'] = hex(struct.unpack('<H', f.read(2))[0])
        f.seek(pe+24+68)
        report['pe_subsystem'] = struct.unpack('<H', f.read(2))[0]
    assert report['pe_machine'] == '0x8664' and report['pe_subsystem'] == 2
    print('Running bundled three-band model verification...', flush=True)
    before = time.perf_counter()
    task = subprocess.Popen([str(exe), '--verify-model', str(output/'model')], **startup)
    try:
        status = task.wait(timeout=300)
    except subprocess.TimeoutExpired:
        subprocess.run(['taskkill', '/PID', str(task.pid), '/T', '/F'], capture_output=True)
        raise
    if status != 0:
        error = output/'model/error.json'
        raise AssertionError(error.read_text(encoding='utf-8') if error.exists() else f'Exe verification exit {status}')
    verification = json.loads((output/'model/verification.json').read_text(encoding='utf-8'))
    assert verification['reference_pairs'] == 1806 and verification['rows'] == 1830
    assert verification['max_abs_prediction_replay_error'] < 1e-10
    report.update(model_verification=verification, full_exe_replay_seconds=time.perf_counter()-before)
    print('Native model replay passed; starting GUI/API checks...', flush=True)
    with socket.socket() as probe:
        probe.bind(('127.0.0.1', 0))
        port = probe.getsockname()[1]
    url = f'http://127.0.0.1:{port}'

    def request(path, body=None):
        if body is None:
            query = urllib.request.Request(url+path)
        else:
            query = urllib.request.Request(url+path, data=json.dumps(body).encode('utf-8'),
                                           headers={'Content-Type': 'application/json'})
        with opener.open(query, timeout=60) as response:
            return response.headers, response.read()

    process = subprocess.Popen([str(exe), '--no-browser', '--port', str(port)], **startup)
    try:
        for _ in range(240):
            if process.poll() is not None:
                raise AssertionError(f'Exe GUI exited with {process.returncode}')
            try:
                health = json.loads(request('/api/health')[1])
                break
            except (OSError, urllib.error.URLError):
                time.sleep(.5)
        else:
            raise AssertionError('Exe GUI did not become ready')
        version = json.loads(request('/api/version')[1])
        assert health['ok'] and version['version'] == '1.1.0' and version['frozen']
        report['version'] = version
        static = ['/', '/app.js', '/styles.css', '/language.js', '/vendor/three.module.js',
                  '/vendor/controls/OrbitControls.js']
        for path in static:
            assert len(request(path)[1]) > 50
        assert 'v1.1' in request('/')[1].decode('utf-8')
        report['static_resources_passed'] = static
        expected = {}
        with (bundle/'all_results.csv').open(encoding='utf-8', newline='') as f:
            for row in csv.DictReader(f):
                if int(row['wavelength_nm']) == 10500:
                    expected[(row['scene'], float(row['vza']), float(row['vaa']))] = row
        manifest = json.loads((bundle/'manifest.json').read_text(encoding='utf-8'))
        matched, hidden, maximum = 0, 0, 0.
        for scene in manifest['scenes']:
            inp = json.loads((bundle/'inputs'/f"{scene['stem']}.json").read_text(encoding='utf-8'))
            config = dict(inp['thermal_config'], spectral_source='direct', wavelengths='10.5',
                          geometry_mode=0, geometry_pair_mode='cross',
                          vza=';'.join(map(str, range(0, 76, 5))), vaa='0;90;180;270')
            for name in ('soil', 'leaf', 'terrain', 'roof', 'wall', 'street'):
                config['emissivity_'+name] = inp['emissivity'][name+'_e']
            data = json.loads(request('/api/run', {'overrides': config})[1])
            assert data['ok'], data
            for row in data['rows']:
                if float(row['vza']) == 0 and float(row['vaa']) != 0:
                    continue
                ref = expected[(scene['stem'], float(row['vza']), float(row['vaa']))]
                if ref['status'] != 'ok':
                    hidden += 1
                    continue
                current = float(row['brightness_temperature_K'])
                difference = abs(current-float(ref['value']))
                maximum = max(maximum, difference)
                assert math.isfinite(current) and difference < 1e-6, (scene['stem'], row, ref)
                assert abs(current-float(row['brightness_temperature_C'])-273.15) < 1e-6
                matched += 1
            print('API passed:', scene['stem'], flush=True)
        # run.py omits the eight back-facing slope views from GUI output;
        # the bundled scientific table retains their explicit status rows.
        assert matched == 602 and hidden == 0
        headers, archive = request('/api/export', {'overrides': config})
        with zipfile.ZipFile(io.BytesIO(archive)) as z:
            assert any(name.endswith('input.csv') for name in z.namelist())
        report.update(api_unique_thermal_predictions=matched, api_backface_rows_returned=hidden,
                      api_hidden_views_filtered=8,
                      api_max_abs_replay_error=maximum, export_zip_passed=True)
        try:
            request('/api/run', [])
        except urllib.error.HTTPError as error:
            assert error.code == 400
        else:
            raise AssertionError('Malformed input was accepted')
        request('/api/shutdown', {})
        assert process.wait(timeout=20) == 0
        report['shutdown_passed'] = True
    finally:
        if process.poll() is None:
            subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'], capture_output=True)
    report['ok'] = True
    (output/'exe_verification.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'model_verification'}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
