# `scripts/` 开发与发布脚本

脚本都应从项目根目录执行：

```powershell
python validation/validate_current_model.py --plots
python scripts/build_installer.py
python scripts/capture_tirt_screenshots.py
```

- `validate_current_model.py`：兼容转发入口，主程序已移至 `validation/validate_current_model.py`；用已保存输入重算 TIRT，并与冻结 HiSTREAM 结果对比，详见 [validation](../validation/README.md)。
- `build_optical_validation_report.py`：从已保存的 2026-10-24 光学配对结果核算误差并生成独立光学报告，不重跑模型。先运行 `--charts`（matplotlib），再运行 `--pdf`（reportlab、Pillow、pypdf）；输出 PDF、Markdown、图件及 CSV。
- `build_installer.py`：构建 Windows 安装包 `release/v<版本>/TiRT-Setup-<版本>-x64.exe`，同时在该目录输出代码压缩包 `TiRT-<版本>-source.zip` 和 PDF 文档 `docs/`。
- `capture_tirt_screenshots.py`：启动本地 GUI，用 Firefox 生成界面截图。

脚本产生的结果放入 `cases/` 或 `docs/`，不要把临时浏览器配置、缓存和构建目录写入源码目录。
