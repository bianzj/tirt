"""Build the TiRT Windows release: installer, source zip and PDF documents.

python scripts/build_installer.py

Result in release/v<version>/:
  TiRT-Setup-<version>-x64.exe   installer (desktop program; also installs the zip and docs)
  TiRT-<version>-source.zip      model code, inputs, data and Markdown documents
  docs/                          user guides, theory, test and project documents (PDF)

Requires PyInstaller and markdown in this Python, Microsoft Edge (PDF printing),
and Node.js with node_modules (npm install).
"""
import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from base.version import __version__  # noqa: E402

BUILD = ROOT / "build" / "installer"
BUNDLE = "validation/reference"

DOCUMENTS = (
    ("README.md", "01_项目说明", "TiRT 项目说明"),
    ("TiRT_user_guide.md", "02_TiRT方向模型操作说明", "TiRT 方向模型操作说明"),
    ("TIRTEB_user_guide.md", "03_TIRTEB能量平衡操作说明", "TIRTEB 能量平衡操作说明"),
    ("TiRT_theory_and_methods.md", "04_TiRT理论与方法", "TiRT 理论与方法"),
    ("TIRTEB_theory_and_methods.md", "05_TIRTEB理论与方法", "TIRTEB 理论与方法"),
    ("TiRT_validation.md", "06_测试文档_HiSTREAM对比", "测试文档：与 HiSTREAM 三维路径追踪对比"),
    ("PROJECT_FILES.md", "07_项目代码说明", "TiRT 项目代码说明"),
)
EDGE = (r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe")

SOURCE_FILES = ("run.py", "tirteb_run.py", "optical_run.py", "input.csv", "input_tirteb.csv", "requirements.txt",
                "environment.yml", "README.md", "PROJECT_FILES.md", "LICENSE", "package.json", "electron-builder.yml",
                "tirt.ico", "start_tirt.bat", "start_tirt.command", ".gitignore")
SOURCE_DIRS = ("base", "gui", "data", "plot", "analysis", "scripts", "docs", "desktop", "validation", "cases")
SOURCE_SKIP = re.compile(r"(__pycache__|\.pyc$|^validation/[^/]+/runs/|^validation/results/|^cases/[^/]+/(?!input\.csv$)|^cases/_)")

CSS = """
@page{size:A4;margin:16mm 15mm}
body{font-family:"Microsoft YaHei","PingFang SC",sans-serif;line-height:1.65;color:#1f2933;font-size:10.5pt}
h1,h2,h3{color:#102a43;line-height:1.3;page-break-after:avoid}h1{border-bottom:2px solid #c62828;padding-bottom:.3em}
a{color:#1565c0}code{background:#f0f4f8;padding:.1em .3em;border-radius:3px;font-size:.9em}
pre{background:#f0f4f8;padding:.7em 1em;white-space:pre-wrap;word-break:break-all;border-radius:4px}pre code{background:none;padding:0}
table{border-collapse:collapse;margin:.8em 0;page-break-inside:avoid}th,td{border:1px solid #bcccdc;padding:.3em .6em}th{background:#e8eef4}
img{max-width:100%;page-break-inside:avoid}
"""


def log(message):
    print(f"[build] {message}", flush=True)


def page(title, body):
    return (f'<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>{title} · TiRT {__version__}</title>'
            f'<style>{CSS}</style></head><body>{body}</body></html>')


def _print_pdf(html, pdf, profile):
    edge = next((path for path in EDGE if Path(path).exists()), None) or shutil.which("msedge")
    if not edge:
        raise SystemExit("Microsoft Edge is required to print the documents to PDF.")
    subprocess.run([edge, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", f"--user-data-dir={profile}",
                    f"--print-to-pdf={pdf}", html.resolve().as_uri()],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=180)
    if not pdf.exists() or pdf.stat().st_size == 0:
        raise SystemExit(f"PDF printing failed: {pdf.name}")


def build_docs(target):
    """User guides, theory, test and project documents as PDF files."""
    import markdown
    from markdown.extensions.toc import slugify_unicode
    target.mkdir(parents=True)
    work = BUILD / "docs_html"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    shutil.copytree(ROOT / "docs/validation/current_model_20261021/figures", work / "validation/current_model_20261021/figures",
                    ignore=shutil.ignore_patterns("*.pdf"))
    pdf_names = {name: f"{stem}.pdf" for name, stem, _ in DOCUMENTS}

    def link(match):
        url = match.group(2)
        name = url.partition("#")[0].split("/")[-1]
        if name in pdf_names:
            url = pdf_names[name]
        elif name == "LICENSE":
            url = "LICENSE.txt"
        return f"{match.group(1)}({url})"

    for name, stem, title in DOCUMENTS:
        source = ROOT / name if name in ("README.md", "PROJECT_FILES.md") else ROOT / "docs" / name
        text = re.sub(r"(\]\s*)\(([^)\s]+)\)", link, source.read_text(encoding="utf-8"))
        html = markdown.markdown(text, extensions=["tables", "fenced_code", "toc"],
                                 extension_configs={"toc": {"slugify": slugify_unicode}})
        page_path = work / f"{stem}.html"
        page_path.write_text(page(title, html), encoding="utf-8")
        _print_pdf(page_path, target / f"{stem}.pdf", work / "edge-profile")
    shutil.copy2(ROOT / "LICENSE", target / "LICENSE.txt")
    shutil.rmtree(work, ignore_errors=True)
    log(f"docs: {len(DOCUMENTS)} PDF files")


def build_source(target):
    target.mkdir(parents=True)
    archive = target / f"TiRT-{__version__}-source.zip"
    prefix = f"TiRT-{__version__}/"
    count = 0
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        paths = [ROOT / name for name in SOURCE_FILES]
        for folder in SOURCE_DIRS:
            paths.extend(sorted(p for p in (ROOT / folder).rglob("*") if p.is_file()))
        for path in paths:
            relative = path.relative_to(ROOT).as_posix()
            if SOURCE_SKIP.search(relative):
                continue
            zf.write(path, prefix + relative)
            count += 1
    log(f"source: {archive.name}, {count} files, {archive.stat().st_size / 1e6:.1f} MB")


def git_commit():
    git = shutil.which("git") or next(iter(sorted(Path(os.environ.get("LOCALAPPDATA", "")).glob(
        "GitHubDesktop/app-*/resources/app/git/cmd/git.exe"), reverse=True)), None)
    if not git:
        return "unknown"
    commit = subprocess.run([str(git), "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run([str(git), "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    return commit + ("+local-changes" if dirty else "")


def build_backend(target):
    import numpy, pandas, scipy
    work = BUILD / "pyinstaller"
    work.mkdir(parents=True, exist_ok=True)
    info = dict(version=__version__, source_commit=git_commit(), python=sys.version.split()[0],
                architecture=platform.machine(), built_at=datetime.now(timezone.utc).isoformat(),
                numeric_versions=dict(numpy=numpy.__version__, scipy=scipy.__version__, pandas=pandas.__version__))
    (work / "build_info.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    parts = tuple(int(x) for x in __version__.split(".")) + (0,)
    (work / "windows_version.txt").write_text(f'''VSVersionInfo(ffi=FixedFileInfo(filevers={parts!r}, prodvers={parts!r},
mask=0x3f, flags=0, OS=0x40004, fileType=0x1, subtype=0, date=(0,0)),
kids=[StringFileInfo([StringTable('040904B0',[
StringStruct('CompanyName','TIRTLab'), StringStruct('FileDescription','TiRT {__version__} service'),
StringStruct('FileVersion','{__version__}'), StringStruct('InternalName','TiRT-service'),
StringStruct('OriginalFilename','TiRT-service.exe'), StringStruct('ProductName','TiRT'),
StringStruct('ProductVersion','{__version__}')])]), VarFileInfo([VarStruct('Translation',[1033,1200])])])''',
                                             encoding="utf-8")
    bundle = ROOT / BUNDLE
    data = [(ROOT / "gui", "gui"), (ROOT / "data", "data"), (ROOT / "input.csv", "."), (ROOT / "input_tirteb.csv", "."),
            (work / "build_info.json", "."), (bundle / "inputs", f"{BUNDLE}/inputs")]
    data += [(bundle / name, BUNDLE) for name in ("all_results.csv", "histream_reference.csv", "accuracy.csv", "manifest.json")]
    # The bundled validation replay records checksums of the model source files.
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    data += [(ROOT / name, Path(name).parent.as_posix()) for name in manifest["model_source_sha256"]]
    command = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--onedir", "--name", "TiRT-service",
               "--icon", str(ROOT / "tirt.ico"), "--version-file", str(work / "windows_version.txt"),
               "--paths", str(ROOT), "--paths", str(ROOT / "validation/streamsim"),
               "--collect-submodules", "base", "--hidden-import", "scripts.validate_current_model",
               "--hidden-import", "run_benchmark", "--distpath", str(target), "--workpath", str(work / "work"),
               "--specpath", str(work)]
    for name in ("matplotlib", "PyQt5", "PyQt6", "PySide6", "tkinter", "IPython"):
        command += ["--exclude-module", name]
    for source, destination in data:
        command += ["--add-data", f"{source}{os.pathsep}{destination}"]
    command.append(str(ROOT / "gui/launcher.py"))
    subprocess.run(command, cwd=ROOT, env=dict(os.environ, MPLBACKEND="Agg"), check=True)
    exe = target / "TiRT-service" / "TiRT-service.exe"
    size = sum(p.stat().st_size for p in exe.parent.rglob("*") if p.is_file())
    log(f"backend: {exe}, folder {size / 1e6:.0f} MB")
    return info


def build_installer(node):
    cli = ROOT / "node_modules/electron-builder/cli.js"
    if not cli.exists():
        raise SystemExit("node_modules/electron-builder is missing; run `npm install` first.")
    subprocess.run([node, str(cli), "--win", "nsis", "--x64", "--publish", "never"], cwd=ROOT, check=True)
    setup = ROOT / "release" / f"v{__version__}" / f"TiRT-Setup-{__version__}-x64.exe"
    log(f"installer: {setup} ({setup.stat().st_size / 1e6:.0f} MB)")
    return setup


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--node", default=shutil.which("node") or "node", help="Node.js executable")
    parser.add_argument("--skip-backend", action="store_true", help="Reuse build/installer/backend from a previous run.")
    args = parser.parse_args()
    if sys.platform != "win32":
        raise SystemExit("Build the Windows installer on Windows.")
    keep = BUILD / "backend" if args.skip_backend else None
    if BUILD.exists():
        for child in BUILD.iterdir():
            if child != keep and child.name != "pyinstaller":
                shutil.rmtree(child) if child.is_dir() else child.unlink()
    BUILD.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "LICENSE", BUILD / "LICENSE.txt")
    build_docs(BUILD / "docs")
    build_source(BUILD / "source")
    if not args.skip_backend:
        build_backend(BUILD / "backend")
    setup = build_installer(args.node)
    # Release folder: installer, source zip and PDF documents side by side.
    shutil.copy2(BUILD / "source" / f"TiRT-{__version__}-source.zip", setup.parent)
    shutil.copytree(BUILD / "docs", setup.parent / "docs", dirs_exist_ok=True)
    for extra in ("win-unpacked", "builder-debug.yml", "latest.yml", "SHA256SUMS.txt", setup.name + ".blockmap"):
        path = setup.parent / extra
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
    log(f"release: {setup.parent}")


if __name__ == "__main__":
    main()
