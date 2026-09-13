# TiRT 项目文件说明

## 项目定位

项目根目录统一维护 TiRT 的辐射传输、能量平衡、时间序列和多角度能力。用户只需要启动一个 TiRT 桌面窗口；窗口内部使用本地 GUI 服务，底层保留批处理兼容入口。

## 目录职责

| 目录 | 职责 | 是否运行必需 |
| --- | --- | --- |
| `base/` | 核心算法、模型、输入解析、能量平衡和兼容示例 | 是 |
| `data/defined/` | 定义参数、太阳/天空光谱、光学系数、气象数据 | 能量平衡功能需要 |
| `data/directions/` | 观测方向输入 | 文件角度模式需要 |
| `data/spectrals/` | TiRT 示例光谱 | 文件光谱模式需要 |
| `cases/` | 案例输入和本地生成结果 | 运行时按需生成 |
| `gui/` | TiRT 统一页面、静态资源和本地服务 | GUI 使用时需要 |
| `desktop/` | Electron 桌面窗口、服务启动和退出清理 | 桌面窗口需要 |
| `plot/` | TiRT 输出绘图 | 绘图时需要 |
| `analysis/` | 实验分析脚本和分析结果 | 可选 |
| `scripts/` | 自动验证、截图和视频工具 | 开发/发布时使用 |
| `tests/` | 单元测试和 GUI 测试 | 开发/发布时使用 |
| `docs/` | 用户手册、理论说明和测试记录 | 发布时建议保留 |

每个主要目录都带有自己的 `README.md`，进入子目录时可以直接看到该目录的边界、入口和产物约定。

## 根目录入口

- `run.py`：TiRT 主入口，读取 `input.csv`。
- `tirteb_run.py`：TiRT 能量平衡时间序列兼容入口，读取 `input_tirteb.csv`。
- `angle_run.py`：TiRT 多角度/多波段兼容入口。
- `input.csv`、`input_tirteb.csv`：两套输入模板，不要互相替换。
- `start_tirt.*`：TiRT 桌面窗口启动脚本；只需启动它。
- `build_tirt_exe.bat`、`electron-builder.yml`：Windows 桌面版打包配置，生成安装包和 `TiRT.exe`。
- `build/TiRT-service.spec`：由构建脚本自动生成的后台服务配置，不需要手工维护。

## 推荐工作流

1. 启动 `start_tirt.bat` 或 `npm run desktop`。
2. 在统一 GUI 中选择方向性模型或能量平衡功能。
3. 批处理时复制 `input.csv` 或 `input_tirteb.csv` 到 `cases/<case_name>/`，按需要修改参数。
4. 保留案例输入和必要的数据文件；图表和 CSV 结果可随时重新生成。

## 输出约定

- 运行结果放在 `cases/<case_name>/` 或调用者指定的目录。
- `output.csv` 保存主结果。
- 能量平衡还会生成 `observations.csv` 和可选的 `intermediate.csv`。
- `cases/` 中的输出文件可以删除并重新生成；Python 缓存、构建目录和发布目录不应提交。

## 验证命令

```bat
cd /d C:\work\tirt
python scripts\verify_run.py
python tests\test_urban_wall_radiation.py
python -m compileall -q .
```

能量平衡回归验证：

```bat
python scripts\verify_tirteb.py --input input_tirteb.csv --output-root cases\_verify_tirteb
```
