"""Desktop launcher: exit codes, error reporting and parent-exit shutdown."""
import io
import json
import subprocess
import sys
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from gui import launcher

ROOT = Path(__file__).resolve().parents[1]


class _FakeCtypes:
    def __init__(self):
        self.calls = []
        self.windll = SimpleNamespace(user32=SimpleNamespace(MessageBoxW=self._box))

    def _box(self, *args):
        self.calls.append(args)
        return 0


def _failing_server(*_args, **_kwargs):
    raise OSError("Ports 8765-8784 are unavailable")


def test_server_only_start_failure_exits_nonzero_without_dialog(monkeypatch, capsys):
    fake = _FakeCtypes()
    monkeypatch.setattr(launcher, "ctypes", fake)
    monkeypatch.setattr(launcher, "_create_server", _failing_server)
    assert launcher.main(["--server-only", "--port", "0"]) == 1
    assert fake.calls == []  # Electron shows the stderr text itself
    assert "could not start" in capsys.readouterr().err


@pytest.mark.skipif(sys.platform != "win32", reason="MessageBox is Windows only")
def test_interactive_start_failure_shows_dialog_and_exits_nonzero(monkeypatch):
    fake = _FakeCtypes()
    monkeypatch.setattr(launcher, "ctypes", fake)
    monkeypatch.setattr(launcher, "_create_server", _failing_server)
    assert launcher.main(["--no-browser", "--port", "0"]) == 1
    assert len(fake.calls) == 1 and "could not start" in fake.calls[0][1]


def test_watch_parent_stops_server_on_stdin_eof():
    stopped = threading.Event()
    server = SimpleNamespace(shutdown=stopped.set)
    launcher._watch_parent(server, io.BytesIO(b""))
    assert stopped.wait(5)


def test_backend_exits_when_parent_closes_stdin(tmp_path):
    """Simulates Electron quitting/crashing: the stdin pipe closes and the backend stops."""
    ready = tmp_path / "ready.json"
    process = subprocess.Popen(
        [sys.executable, str(ROOT / "gui" / "launcher.py"), "--host", "127.0.0.1", "--port", "0",
         "--server-only", "--ready-file", str(ready)],
        cwd=ROOT, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    try:
        deadline = time.time() + 60
        while not ready.exists() and process.poll() is None and time.time() < deadline:
            time.sleep(0.1)
        assert ready.exists(), process.stderr.read().decode("utf-8", "replace") if process.poll() is not None else "no ready file"
        assert json.loads(ready.read_text(encoding="utf-8"))["port"] > 0
        process.stdin.close()
        assert process.wait(timeout=20) == 0
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
