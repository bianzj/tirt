# `scripts/` 开发与发布脚本

脚本都应从项目根目录执行：

```powershell
python scripts/validate_current_model.py --plots
python scripts/build_tirt_release.py
python scripts/capture_tirt_screenshots.py
```

- `validate_current_model.py`：用 `docs/validation/current_model_20261021/` 中的输入重算与 HiSTREAM 的对比，并校验与包内结果一致。
- `build_tirt_release.py`：构建 Windows 单文件 `TiRT.exe`。
- `capture_tirt_screenshots.py`：启动本地 GUI，用 Firefox 生成界面截图。

脚本产生的结果放入 `cases/` 或 `docs/`，不要把临时浏览器配置、缓存和构建目录写入源码目录。
