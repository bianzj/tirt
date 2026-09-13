# `scripts/` 开发与发布脚本

脚本都应从项目根目录执行：

```powershell
python scripts/verify_run.py
python scripts/verify_tirteb.py --input input_tirteb.csv --output-root cases/_verify_tirteb
```

- `verify_*.py`：回归验证和快速烟囱测试。
- `capture_*.py`：启动本地 GUI 并生成截图。
- `make_*.py`：根据截图和结果生成演示视频。

脚本产生的结果放入 `cases/` 或 `docs/`，不要把临时浏览器配置、缓存和构建目录写入源码目录。
