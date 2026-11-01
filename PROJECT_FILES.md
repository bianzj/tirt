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
| `cases/` | 案例输入和本地生成结果（CSV 运行默认写入 `cases/<case_name>/`） | 运行时按需生成 |
| `gui/` | TiRT 统一页面、静态资源和本地服务 | GUI 使用时需要 |
| `desktop/` | Electron 桌面窗口、服务启动和退出清理 | 桌面窗口需要 |
| `plot/` | TiRT 输出绘图 | 绘图时需要 |
| `analysis/` | 实验分析脚本和分析结果 | 可选 |
| `scripts/` | HiSTREAM 验证复算、打包和截图工具 | 开发/发布时使用 |
| `docs/` | 操作说明、理论说明、测试文档（HiSTREAM 对比）与验证包 | 发布时建议保留 |

每个主要目录都带有自己的 `README.md`，进入子目录时可以直接看到该目录的边界、入口和产物约定。

## 根目录入口

- `optical_run.py`：TiRT 光学运行入口，用 `--input <场景>.json` 指定案例；输出多光谱/多角度 CSV 和两张图。
- `thermal_run.py`：读取热红外 JSON，默认多 LAI、多光谱和多角度 CSV 与两张图；直接指定组分温度，不求能量平衡，不带 `--input` 时读取 `input.csv`。
- `run.py`：方向性兼容入口与共用后端，同一张 `output.csv` 可同时包含短波 `reflectance_brf` 和 10.5 µm 亮温，已有调用保持可用。
- `base/scenario_runner.py`：多 LAI、光谱/角度扫描、结果导出与可选绘图；不依赖 GUI 或案例脚本。
- `tirteb_run.py`：能量平衡入口，默认温度/通量时间序列 CSV 与两张图；不带 `--input` 时读取 `input_tirteb.csv`。
- `base/energy_runner.py`：能量平衡 JSON 转换、引用路径解析、时间序列导出和绘图；数值求解仍使用原后端。
- `base/single_point.py`：单场景能量平衡（一个方向、一个波段），`python -m base.single_point`。
- `input.csv`、`input_tirteb.csv`：两个输入模板，使用同一套参数名称和含义（共有参数由 `base/inputio.py` 解析）。`input.csv` 额外给出组分温度，`input_tirteb.csv` 额外给出气象、能量平衡和短波光学参数。
- `start_tirt.*`：TiRT 桌面窗口启动脚本；只需启动它。
- `scripts/build_installer.py`、`electron-builder.yml`、`desktop/installer.nsh`：Windows 安装包构建；中间文件写入 `build/installer/`，安装包写入 `release/v<版本>/`。

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

## 验证

项目不再附带 `examples/` 场景示例；模型复算与精度检查统一使用 `validation/`。
与 HiSTREAM 的对比见[测试文档](docs/TiRT_validation.md)；多山脊 DEM 的对照见 [validation/ridge_20261029](validation/ridge_20261029/README.md)。复算命令：

```bat
cd /d C:\work\tirt
python validation\validate_current_model.py --plots
```
