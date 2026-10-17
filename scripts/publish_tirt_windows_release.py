"""Publish verified TiRT Windows assets using the existing Git credential helper."""
import argparse
import hashlib
import json
import mimetypes
import os
import subprocess
import time
from pathlib import Path
from urllib.parse import urlsplit

import requests


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


class ProgressFile:
    def __init__(self, path):
        self.path = path
        self.file = path.open('rb')
        self.size = path.stat().st_size
        self.last = time.monotonic()

    def __len__(self):
        return self.size

    def tell(self):
        return self.file.tell()

    def seek(self, *args):
        return self.file.seek(*args)

    def read(self, size=-1):
        block = self.file.read(size)
        if time.monotonic()-self.last > 20:
            print(f'{self.path.name}: {self.tell()/1024**2:.1f}/{self.size/1024**2:.1f} MiB', flush=True)
            self.last = time.monotonic()
        return block


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--assets', type=Path, required=True)
    args = parser.parse_args()
    repository, folder = args.repository.resolve(), args.assets.resolve()
    info = json.loads((folder/'build_info.json').read_text(encoding='utf-8'))
    assert info['version'] == '1.1.2'
    verification = json.loads((folder/'exe_verification.json').read_text(encoding='utf-8'))
    assert verification['dynamic_verification']['uploads_passed']
    assert len(verification['dynamic_verification']['cases']) == 6
    assert verification['ok'] and verification['model_verification']['reference_pairs'] == 1806
    assert sha(folder/'TiRT.exe') == info['sha256'] == verification['sha256']
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repository, text=True).strip()
    built_from = info['source_commit']
    subprocess.run(['git', 'merge-base', '--is-ancestor', built_from, head], cwd=repository, check=True)
    changes = subprocess.check_output(['git', 'diff', '--name-only', built_from, head],
                                      cwd=repository, text=True).splitlines()
    assert set(changes) <= {'scripts/verify_tirt_exe_release.py'}, 'Packaged program source has changed.'
    env = dict(os.environ, GCM_INTERACTIVE='Never', GIT_TERMINAL_PROMPT='0')
    credential = subprocess.run(['git', 'credential', 'fill'], cwd=repository,
        input='protocol=https\nhost=github.com\n\n', text=True, capture_output=True, check=True, env=env)
    values = dict(line.split('=', 1) for line in credential.stdout.splitlines() if '=' in line)
    token = values.get('password')
    assert token, 'Existing GitHub credential is unavailable.'
    session = requests.Session()
    session.proxies.update({'https': 'http://127.0.0.1:7897', 'http': 'http://127.0.0.1:7897'})
    session.headers.update({'Authorization': 'Bearer '+token, 'Accept': 'application/vnd.github+json',
                            'X-GitHub-Api-Version': '2022-11-28'})
    api = 'https://api.github.com/repos/bianzj/tirt'
    response = session.get(api+'/releases/tags/v1.1.2', timeout=30)
    if response.status_code == 404:
        notes = (folder/'RELEASE_NOTES.md').read_text(encoding='utf-8')
        response = session.post(api+'/releases', json=dict(tag_name='v1.1.2', target_commitish=head,
            name='TiRT 1.1.2', body=notes, draft=True, prerelease=False), timeout=30)
        response.raise_for_status()
        release = response.json()
    else:
        response.raise_for_status()
        release = response.json()
        if not release['draft']:
            raise SystemExit('v1.1.2 is already published; refusing to replace a published release.')
    receipt = folder/'github_release_receipt.json'
    receipt.write_text(json.dumps({'release_id': release['id'], 'draft': True,
        'tag': 'v1.1.2', 'source_commit': built_from, 'release_commit': head}, indent=2), encoding='utf-8')
    upload = release['upload_url'].split('{', 1)[0]
    assert urlsplit(upload).hostname == 'uploads.github.com'
    uploaded = {a['name']: a for a in release['assets']}
    assets = [folder/name for name in ('TiRT.exe', 'build_info.json', 'exe_verification.json',
                                       'TiRT-1.1.2-verification.zip', 'SHA256SUMS.txt')]
    verified = []
    for path in assets:
        expected = sha(path)
        existing = uploaded.get(path.name)
        if existing and existing.get('state') == 'uploaded':
            assert existing['size'] == path.stat().st_size
            assert existing.get('digest') in (None, 'sha256:'+expected)
            asset = existing
        else:
            print('Uploading', path.name, path.stat().st_size, 'bytes', flush=True)
            stream = ProgressFile(path)
            try:
                response = session.post(upload, params={'name': path.name}, data=stream,
                    headers={'Content-Type': mimetypes.guess_type(path.name)[0] or 'application/octet-stream',
                             'Content-Length': str(path.stat().st_size)}, timeout=(30, 1200))
            finally:
                stream.file.close()
            response.raise_for_status()
            asset = response.json()
            assert asset['state'] == 'uploaded' and asset['size'] == path.stat().st_size
            assert asset.get('digest') in (None, 'sha256:'+expected)
        verified.append(dict(name=path.name, size=asset['size'], sha256=expected,
                              remote_digest=asset.get('digest'), url=asset['browser_download_url']))
        print('Verified uploaded asset:', path.name, flush=True)
    response = session.patch(api+f"/releases/{release['id']}",
                             json={'draft': False, 'prerelease': False, 'make_latest': 'true'}, timeout=30)
    response.raise_for_status()
    result = response.json()
    assert not result['draft'] and result['tag_name'] == 'v1.1.2'
    latest = session.get(api+'/releases/latest', timeout=30)
    latest.raise_for_status()
    assert latest.json()['id'] == result['id']
    receipt.write_text(json.dumps(dict(release_id=result['id'], release_url=result['html_url'],
        tag='v1.1.2', source_commit=built_from, release_commit=head, published=True, latest=True, assets=verified),
        ensure_ascii=False, indent=2), encoding='utf-8')
    print('PUBLISHED', result['html_url'], flush=True)


if __name__ == '__main__':
    main()
