"""Build a versioned, standalone Windows TiRT executable from this checkout."""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from base.version import __version__


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dist', type=Path, default=ROOT / 'release/v1.1')
    args = parser.parse_args()
    if sys.platform != 'win32':
        raise SystemExit('Build the Windows executable on Windows.')
    dist = args.dist.resolve()
    dist.mkdir(parents=True, exist_ok=True)
    metadata = ROOT / 'build/tirt_release_metadata'
    metadata.mkdir(parents=True, exist_ok=True)
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    info = dict(version=__version__, source_commit=commit,
                python=sys.version.split()[0], architecture=platform.machine(),
                built_at=datetime.now(timezone.utc).isoformat())
    import numpy
    import scipy
    import pandas
    info['numeric_versions'] = dict(numpy=numpy.__version__, scipy=scipy.__version__, pandas=pandas.__version__)
    python_root = Path(sys.executable).parent
    libdir = python_root/'Library/bin'
    # oneMKL loads CPU/threading kernels at runtime; static PE dependency scans
    # otherwise collect only mkl_rt, producing an exe that fails on clean PCs.
    binaries = [p for p in libdir.glob('mkl*.dll')
                if not any(name in p.name for name in ('blacs', 'scalapack', 'cdft'))]
    for name in ('libiomp5md.dll', 'ffi.dll'):
        if (libdir/name).exists():
            binaries.append(libdir/name)
    info['extra_runtime_libraries'] = [p.name for p in binaries]
    (metadata / 'build_info.json').write_text(json.dumps(info, indent=2), encoding='utf-8')
    parts = tuple(int(x) for x in __version__.split('.')) + (0,)
    version = f'''VSVersionInfo(ffi=FixedFileInfo(filevers={parts!r}, prodvers={parts!r},
mask=0x3f, flags=0, OS=0x40004, fileType=0x1, subtype=0, date=(0,0)),
kids=[StringFileInfo([StringTable('040904B0',[
StringStruct('CompanyName','TIRTLab'), StringStruct('FileDescription','TiRT 1.1'),
StringStruct('FileVersion','{__version__}'), StringStruct('InternalName','TiRT'),
StringStruct('OriginalFilename','TiRT.exe'), StringStruct('ProductName','TiRT'),
StringStruct('ProductVersion','{__version__}')])]), VarFileInfo([VarStruct('Translation',[1033,1200])])])'''
    (metadata / 'windows_version.txt').write_text(version, encoding='utf-8')
    env = dict(os.environ, MPLBACKEND='Agg')
    env['PATH'] = os.pathsep.join([str(python_root), str(python_root/'Library/bin'),
                                  str(python_root/'Scripts'), env.get('PATH', '')])
    bundle = ROOT / 'docs/validation/current_model_20261015'
    data = [(str(ROOT/'gui'), 'gui'), (str(ROOT/'data'), 'data'), (str(ROOT/'input.csv'), '.'),
            (str(metadata/'build_info.json'), '.'),
            (str(bundle/'inputs'), 'docs/validation/current_model_20261015/inputs')]
    for name in ('all_results.csv', 'histream_reference.csv', 'accuracy.csv', 'manifest.json'):
        data.append((str(bundle/name), 'docs/validation/current_model_20261015'))
    # The numerical replay records checksums of the actual model source files.
    # Include these alongside the PYZ modules so that check also works frozen.
    model_manifest = json.loads((bundle/'manifest.json').read_text(encoding='utf-8'))
    for name in model_manifest['model_source_sha256']:
        data.append((str(ROOT/name), Path(name).parent.as_posix()))
    command = [sys.executable, '-m', 'PyInstaller', '--noconfirm', '--clean', '--onefile', '--windowed',
               '--name', 'TiRT', '--icon', str(ROOT/'tirt.ico'), '--version-file', str(metadata/'windows_version.txt'),
               '--paths', str(ROOT), '--paths', str(ROOT/'benchmark/streamsim'),
               '--hidden-import', 'scripts.validate_current_model', '--hidden-import', 'run_benchmark',
               '--distpath', str(dist), '--workpath', str(ROOT/'build/tirt_release'),
               '--specpath', str(metadata)]
    for name in ('matplotlib', 'PyQt6', 'PySide6', 'tkinter'):
        command.extend(['--exclude-module', name])
    for source, target in data:
        command.extend(['--add-data', f'{source}{os.pathsep}{target}'])
    for source in binaries:
        command.extend(['--add-binary', f'{source}{os.pathsep}.'])
    command.append(str(ROOT/'gui/launcher.py'))
    subprocess.run(command, cwd=ROOT, env=env, check=True)
    exe = dist/'TiRT.exe'
    info.update(executable=exe.name, size=exe.stat().st_size,
                sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
    (dist/'build_info.json').write_text(json.dumps(info, indent=2), encoding='utf-8')
    print(json.dumps(info, indent=2))


if __name__ == '__main__':
    main()
