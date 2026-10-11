# `scripts/` 开发与发布脚本

脚本都应从项目根目录执行：

```powershell
python scripts/verify_run.py
python scripts/verify_tirteb.py --input input_tirteb.csv --output-root cases/_verify_tirteb
python scripts/capture_tirt_screenshots.py
```

- `verify_*.py`：回归验证和快速冒烟测试。
- `capture_tirt_screenshots.py`：启动本地 GUI，用 Firefox 生成 `docs/screenshots/` 中的 README 截图。

脚本产生的结果放入 `cases/` 或 `docs/`，不要把临时浏览器配置、缓存和构建目录写入源码目录。
