# `tests/` 测试

测试覆盖核心物理、TiRT GUI、能量平衡 GUI 和页面端到端流程。

```powershell
python tests/test_urban_wall_radiation.py
python tests/test_tirt_gui.py --mode source --port 8765
python tests/test_tirteb_gui.py --mode source --port 8901
```

完整回归优先使用 `scripts/verify_run.py` 和 `scripts/verify_tirteb.py`。页面测试需要 Playwright 浏览器；没有浏览器时仍可运行核心和 API 测试。
