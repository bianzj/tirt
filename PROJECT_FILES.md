# TiRT 文件整理

## 必须保留

- `run.py`：统一计算入口，提供 `run(input_path, output_path, **overrides)` 和命令行运行。
- `input.csv`：默认输入模板；GUI 和命令行都依赖它。
- `base/`：核心模型、物理函数、输入解析和工具函数。
- `data/`：方向文件与光谱示例文件；`input.csv` 的文件模式会引用这里。
- `gui/`：浏览器 GUI、静态资源、本地服务和 exe launcher。
- `requirements.txt`、`environment.yml`：Python/Conda 依赖说明。
- `start_tirt.bat`、`start_tirt.command`：不打包时的一键启动脚本。
- `build_tirt_exe.bat`、`TiRT.spec`：Windows GUI exe 构建入口与 PyInstaller 配置。
- `scripts/verify_run.py`：验证 `run.py`、`input.csv`、模型组合和文件输入模式。

## 可选保留

- `plot/`：读取 `run.py` 输出并生成极坐标/主平面图；如果只使用 GUI，可选。
- `base/run_tirt_*.py`：旧示例入口；主流程已有 `run.py`，这些文件适合作为示例和回归参考。
- `.vscode/`、`.idea/`：本机 IDE 配置；协作发布包中可不带。
- `dist/TiRT.exe`：单文件发布版；已内置默认 `gui/`、`data/`、`input.csv`。

## 可删或不应提交

- `build/`：PyInstaller 中间构建目录，可重新生成。
- `cases/`：运行输出目录，可重新生成。
- `__pycache__/`、`base/__pycache__/`、`gui/__pycache__/`、`plot/__pycache__/`、`rt/__pycache__/`：Python 缓存。
- `rt/`：当前只包含 `__pycache__`，没有源码文件，可删除。
- `.DS_Store`：macOS 缓存文件。
- `.Rhistory`：R 会话历史；当前已被 Git 跟踪，后续建议从版本库移除。

## 验证命令

```bat
cd C:\work\tirt
C:\work\miniconda\python.exe scripts\verify_run.py
C:\work\miniconda\python.exe run.py --output cases\_verify\default_cli\output.csv
C:\work\miniconda\python.exe -m plot.run_plot --input cases\_verify\default_cli\output.csv --mode all --wavelength 10.5 --output-dir cases\_verify\default_cli\plots
```
