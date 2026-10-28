# TiRT

TiRT (Terrestrial Inhomogeneous Radiative Transfer) is a land-surface radiative transfer model for vegetation, terrain, urban surfaces, and coupled scenes.

文档：

- 项目代码说明：[PROJECT_FILES.md](PROJECT_FILES.md)
- 操作说明：[TiRT 操作说明](docs/TiRT_user_guide.md)、[TIRTEB 能量平衡操作说明](docs/TIRTEB_user_guide.md)
- 理论说明：[TiRT 理论与方法](docs/TiRT_theory_and_methods.md)、[TIRTEB 理论与方法](docs/TIRTEB_theory_and_methods.md)
- 验证报告：[GitHub 主页热红外验证](docs/TiRT_validation.md) · [PDF](docs/TIRT_热红外验证报告.pdf)
- 光学验证：[10 场景红光与近红外验证](docs/TIRT_光学验证报告.md) · [PDF](docs/TIRT_光学验证报告.pdf)，研究来源：闻建光团队。

## TiRT 1.2.0 Windows 安装版

安装包 `TiRT-Setup-1.2.0-x64.exe` 安装后得到 TiRT 桌面程序（无需安装 Python），并在安装目录附带：

- `docs/`：操作说明、理论说明、测试文档和项目代码说明（PDF，开始菜单有 `TiRT 文档` 快捷方式）；
- `TiRT-1.2.0-source.zip`：代码压缩包（模型源代码、输入模板、数据与 Markdown 文档）。

发布目录 `release/v1.2.0/` 中同样是这三项：安装程序、代码压缩包和 `docs/`。

1.2.0 采用参数化标准（树冠、垄行闭式四组分比例），TIRTEB 的热红外观测直接调用同一 TiRT 模型。

```powershell
# 从源码构建安装包（需要 PyInstaller、Node.js 与 node_modules 中的 electron-builder）
python scripts/build_installer.py
# 使用安装后的后台服务复算内置模型验证
& "$env:LOCALAPPDATA\Programs\TiRT\resources\backend\TiRT-service.exe" --verify-model .\verification
```

## GitHub 主页热红外验证（2026-10-11）

热红外验证报告整理自 GitHub `master` 提交 `86c3829` 的主页：保留 11 个场景的原始对比图，以及 TIRT、一维均质基准相对 HiSTREAM 的误差统计，共 696 条统计记录。热红外研究来源为杜永明团队。见[报告正文](docs/TiRT_validation.md)与 [PDF](docs/TIRT_热红外验证报告.pdf)。本报告整理已发布结果，未重新运行模型。

当前代码的独立复算入口为 `python validation/validate_current_model.py`，加 `--plots` 可绘图；说明见 [validation](validation/README.md)。其结果与上述历史主页报告分别记录。

## 一、简介

TiRT（Terrestrial Inhomogeneous Radiative Transfer，陆表非均质辐射传输）用于陆表非均质辐射传输建模，光学与热红外共用同一套地表、植被和城市结构模型。热红外链计算不同观测角度下的方向亮温、辐射亮度和有效发射率，光学链计算方向反射因子（BRF）。

TiRT 包含三条计算链：

- **光学方向性辐射传输（`optical_run.py`）**：读取 `--input` 指定的场景 JSON，计算不同 LAI 下的多波段和多角度反射因子，保存数据与两张图。
- **热红外方向性辐射传输（`thermal_run.py`）**：读取 JSON 或 CSV，给定组分温度，计算多 LAI、多波段、多角度亮温和辐亮度，默认保存两张图与 CSV。`run.py` 保留为兼容入口与共用后端。
- **能量平衡时间序列（TIRTEB，`tirteb_run.py`）**：由气象数据驱动，逐时间步求解土壤、叶片和城市组分的温度与 Rn、LE、H、G，再用同一套 TiRT 辐射传输计算方向亮温（直接调用 `run.py` 的模型与参数化标准，同一场景、相同组分温度下两者辐亮度一致；城市植被由所选冠层模型（均质、垄行或离散树冠）直接计算，不另设覆盖比例）。两者共用一套参数名称和含义（`input.csv` 与 `input_tirteb.csv`）。

方向性入口默认使用 **0.67 µm（红光）、0.86 µm（近红外）、10.5 µm（热红外）** 三个波段：前两个短波波段输出反射率，10.5 µm 输出亮温。短波默认按 80% 直射 + 20% 漫射（`shortwave_diffuse_fraction=0.2`）；设为 0 可回到纯直射 BRF（验证报告使用的口径）。静态方向页只保留场景和波段两个选择，选到短波显示反射率，选到 10.5 µm 显示亮温 °C。

LAI 输入支持单个值或分号分隔的列表（如 `1;3;5`）。多 LAI 时，主平面曲线叠加各 LAI，极坐标图和指标用 LAI 下拉切换，`observations.csv` 增加 `lai` 列。

坡面裸地/均质/垄行的短波 BRF 采用局地平面变换；复合地形的非森林组合使用复合地形组分比例计算。复合地形森林的光学算法取决于 `canopy_closure`：默认 `standard`（以及 `upgrade`）使用标准复合地形闭合（`base/optical_terrain_standard.py`），此时 `optical_terrain_method` 不起作用（`effective_slab` 除外）；`canopy_closure=legacy` 时默认 `reciprocal_exchange`（项目内的 4 状态互易闭合），可用 `optical_terrain_method=facet_exchange` 切换为显式面元交换（慢约两个数量级，作为慢速几何参照）。当前精度见[光学验证报告](docs/TIRT_光学验证报告.md)。

### 适用场景

- 均质植被
- 垄行作物
- 离散树冠
- 坡面植被
- 复杂地形
- 城市建筑
- 城市建筑与街道植被耦合场景

### 核心模型

- `Hom`: 均质植被解析模型
- `Row`: 垄行作物解析模型
- `Crown`: 离散树冠解析模型
- `Slope_Veg`: 坡面植被模型
- `Terrain`: 地形模型
- `TerrainForest`: 地形-植被耦合模型
- `Urban`: 城市建筑模型
- `Urban_Veg`: 城市-植被耦合模型

能量平衡部分的方法（组分通量、冠层辐射、空气动力学和气孔阻抗）见 [TIRTEB 理论与方法说明](docs/TIRTEB_theory_and_methods.md)。

### 环境安装

建议使用 Python 3.10+。

```powershell
cd C:\work\tirt
py -m venv .venv
.\.venv\Scripts\activate
pip install -r requirements.txt
```

依赖只有 numpy、scipy、matplotlib、pandas 和 pillow（后者仅 HiSTREAM 对比脚本使用）。

### 统一运行入口

光学入口必须用 `--input` 指定场景 JSON；热红外和能量平衡入口默认读取根目录的
`input.csv`、`input_tirteb.csv`，CSV 结果写入 `cases/<case_name>/`。CSV 导出与两张图由 Python 控制：

```bash
python thermal_run.py
python tirteb_run.py
```

也可以在 Python 中调用外层函数，或通过命令行切换模型：

```python
from thermal_run import run

output_path = run("input.csv", "cases/hom_plane/output.csv")
print(output_path)
```

CSV 输入返回 `output.csv` 路径，不绘图；`output_dir`、`plots`、`show` 只用于 JSON 场景（如 `run("scene.json", output_dir="results/scene", plots=False)`）。

```bash
python thermal_run.py --input input.csv --output cases/row_demo/output.csv --model row
```

统一入口支持植被与地表的组合模型：`plane`、`slope`、`terrain`、`urban` 地表可分别组合 `bare`、`hom`、`row`、`crown` 植被。新案例参数放在 JSON；原 `input.csv` 和 CSV 输出接口保持兼容。

能量平衡时间序列：

```bash
python tirteb_run.py --input input_tirteb.csv --output cases/my_case/output.csv
```

能量平衡 JSON 输入（`--input <scene.json>`）默认输出温度和能量通量的时间序列图与对应 CSV。
`plots=False` 只关闭绘图，仍导出 CSV。

| 需求 | 入口 | 输入 | 输出 |
| --- | --- | --- | --- |
| TiRT 光学方向反射率 | `optical_run.py --input <scene.json>` | 光学 JSON 配置 | `spectrum.csv`、`angles.csv`、两张 PNG |
| TiRT 热红外方向亮温、发射率 | `thermal_run.py` | 热红外 JSON（兼容 CSV） | `spectrum.csv`、`angles.csv`、两张 PNG；CSV 模式保留原输出 |
| TiRT 桌面界面 | `desktop/main.mjs` | `input.csv`、`input_tirteb.csv` | Electron 窗口 |
| TiRT 能量平衡批处理 | `tirteb_run.py` | 能量平衡 JSON（兼容 CSV） | `output.csv`、`observations.csv`；JSON 模式增加温度/通量 CSV 和两张 PNG |
| 单场景能量平衡（一个方向、一个波段） | `python -m base.single_point --vza 0 --vaa 0 --wavelength 10.5` | `input_tirteb.csv` | `cases/single_point/output.csv`（逐时间步组分温度、通量与亮温） |

### input.csv 结构

输入分组采用 `section,key,value,description` 格式：

- `[surface]`：`surface_model` 选择平面、单一坡、复合地形或城区地表，`vegetation_model` 选择裸地、均质植被、垄行或森林。
- `[vegetation_structure]`、`[row_structure]`、`[crown_structure]`：分别配置 LAI、热点、垄行和树冠结构。
- `[slope]`：配置单一坡的坡度和坡向；`[terrain]` 配置复合地形层，可用分号提供多个高度、半径和密度参数。
- `[optical]`：`optical_terrain_method` 仅在 `canopy_closure=legacy` 时控制复合地形森林的光学算法，默认 `reciprocal_exchange`（项目内的 4 状态互易闭合，不是外部命名模型），可改为 `facet_exchange`（显式面元辐射度，慢速参照）；`standard`/`upgrade` 使用标准复合地形闭合。热红外 `forest_transport=reciprocal_exchange` 同样要求 `canopy_closure=legacy`。
- `[spectral]`：光学和热红外共用一张波段参数表，按行依次为 `wavelengths`（波段，µm）、`leaf_reflectance`（叶片反射率）、`leaf_transmittance`（叶片透过率）、`soil_reflectance`（土壤反射率）、`street_reflectance`（街道反射率）、`wall_reflectance`（墙壁反射率）、`roof_reflectance`（屋顶反射率），每行用分号给出每个波段一个值。小于 2.5 µm 的波段直接用这些反射率/透过率计算 BRF；大于等于 2.5 µm 的波段按基尔霍夫定律换算发射率：叶片 `ε = 1 − ρ − τ`，土壤（含复合地形山体）、街道、墙壁、屋顶 `ε = 1 − ρ`。输出 CSV 的 `emissivity_*` 列即换算结果，可直接核对。`spectral_source=file` 时从 `spectrals_file`（默认 [optical_properties.txt](data/spectrals/optical_properties.txt)）读取同样 7 行的表，每行前面可以带文字标签。场景用不到的行可以留空（例如非城市场景的街道/墙壁/屋顶）。
- `[soil_anisotropy]`：选择 `bare` 时，`hapke_enable=1` 启用 Hapke/Thapke 土壤方向性修正（以该波段土壤发射率为半球发射率），`hapke_K` 控制各向异性强度。
- `[thermal]`：配置土壤、叶片、地形、屋顶、墙壁和街道的日照/阴影温度。`sky_temperature_k` 可选：留空时输出仅含发射（黑天空，与验证一致）；填写天空等效亮温后，按 `(1 − 有效发射率) × B(T_sky)` 加入天空下行辐射的反射。
- `[row_structure]`、`[crown_structure]` 中的 `thermal_row_geometry`、`thermal_crown_geometry`、`forest_transport` 默认即参数化标准（`parametric`、`parametric`、`legacy`），一般无需修改；树冠参数化需要冠底高度 `crown_height`。`crown_rad_a` 为水平半径、`crown_rad_b` 为垂直半轴；树冠的 `hspot` 取叶片大小 / 树冠高度。
- `[geometry]`：通过 `geometry_mode` 选择观测角度来源：`0` 为 `vza/vaa` 给定角度，`1` 自动生成主平面，`2` 从 TXT 文件读取，`3` 以一个天顶角步长同时生成太阳主平面和垂直太阳主平面，并追加文件或手动的一一对应观测角度。

能量平衡模板 `input_tirteb.csv` 使用同一套参数名称和含义（地表、坡面、地形、建筑、垄行、树冠、波段参数表、观测几何），其中波段参数表的格可以填 `-`，表示短波取 PROSPECT/BSM 或实测光谱、热红外取默认值；能量平衡的长波收支使用表中热红外波段的发射率。组分温度由能量平衡求出，热红外方向辐射直接调用 TiRT 的辐射传输，见 [TIRTEB 理论与方法说明](docs/TIRTEB_theory_and_methods.md#参数体系)。

几何文件每行支持 `vza vaa`，也支持 `vza vaa sza saa`；角度单位为度，注释行以 `#` 开头。模式 1 使用相对太阳方位角 `0/90/180/270`，可通过 `principal_vza` 或 `principal_vza_step`、`principal_vza_max` 设置天顶角；设置 `geometry_view=parallel` 或 `geometry_view=perpendicular` 可只生成对应的两个方向。

项目提供了半球观测方向文件 `data/directions/hemisphere_directions.txt`：方位角间隔 30°、天顶角间隔 10°，覆盖 0° 到 90° 天顶角。使用该文件时设置：

```text
geometry,geometry_mode,2,read TXT directions
geometry,geometry_file,data/directions/hemisphere_directions.txt,hemisphere directions
```

### 方向性绘图

`plot/` 用于读取统一入口生成的 `output.csv`，绘制热红外方向性结果。默认绘制摄氏度亮温，也可以选择开尔文亮温、辐亮度或其他数值列：

```bash
python -m plot.run_plot \
    --input cases/file_geometry_spectrum/output.csv \
    --mode all \
    --wavelength 10.5
```

可选模式为：

- `polar`：极坐标方向图，角度为观测方位角 `vaa`，半径为观测天顶角 `vza`。
- `parallel`：太阳主平面图，使用相对太阳方位角 `raa=0/180` 的方向。
- `perpendicular`：垂直太阳主平面图，使用相对太阳方位角 `raa=90/270` 的方向。
- `all`：一次输出以上三类图；如果输入几何没有对应方向，则跳过该类图。

默认图像写入 `cases/<case_name>/plots/`，也可以指定目录和绘制量：

```bash
python -m plot.run_plot \
    --input cases/file_geometry_spectrum/output.csv \
    --mode polar \
    --wavelength 10.5 \
    --quantity brightness_temperature_K \
    --output-dir cases/file_geometry_spectrum/plots_K
```

绘图依赖 `matplotlib`，安装项目依赖后即可使用。

### 基本用法

```python
import numpy as np
from base.hom import Hom

lai = 0.5
hspot = 0.15
wavelength = 10.5
emissivity_leaf = 0.985
emissivity_soil = 0.955

temperature_leaf_sunlit = 303
temperature_leaf_shaded = 300
temperature_soil_sunlit = 320
temperature_soil_shaded = 305

sza = 25
vza = np.hstack([np.linspace(50, 1, 50), np.linspace(0, 50, 51)])
raa = np.zeros_like(vza)

hom = Hom()
hom.set_structure(lai, hspot)
hom.set_optical(wavelength, emissivity_soil, emissivity_leaf)
hom.set_thermal(
    temperature_soil_sunlit,
    temperature_soil_shaded,
    temperature_leaf_sunlit,
    temperature_leaf_shaded,
)
hom.set_angle(vza, sza, raa)

brightness_temperature = hom.run()
```

### 目录结构

```text
tirt/
├── optical_run.py          # TiRT 光学运行入口
├── thermal_run.py          # TiRT 热红外运行入口
├── run.py                  # 热红外兼容入口与共用后端
├── tirteb_run.py           # TiRT 能量平衡批处理入口
├── input.csv               # TiRT 默认输入模板
├── input_tirteb.csv        # 能量平衡默认输入模板
├── base/                   # 统一 Python 核心包（保持 base 导入兼容）
│   ├── gap.py              # 间隙概率与几何路径
│   ├── hotspot.py          # 热点效应
│   ├── proportion.py       # 四组分方向比例
│   ├── scatter.py          # VNIR/TIR 散射
│   ├── hom/row/crown*.py  # 均质、垄行、树冠及体素模型
│   ├── terrain*.py         # 坡面、复合地形及地形植被模型
│   ├── urban*.py           # 城市建筑及城市植被模型
│   ├── inputio.py          # TiRT 与 TIRTEB 共用的参数解析
│   ├── radiation.py        # 统一辐射准备、合并和诊断
│   ├── rt_transfer.py      # 能量平衡的冠层短波、长波传输
│   ├── aerodynamic.py     # 空气动力学阻抗
│   ├── biochemical.py      # 光合/气孔阻抗
│   ├── evapo.py/eb.py      # 蒸散、能量平衡和热量更新
│   ├── io.py/modelio.py    # 能量平衡输入、光谱和结构解析
│   ├── structs.py          # 能量平衡数据结构
│   └── single_point.py     # 单场景能量平衡入口
├── data/
│   ├── defined/            # 光谱系数、太阳光谱和气象数据
│   ├── directions/         # 观测方向文件
│   └── spectrals/          # TiRT 示例观测光谱
├── validation/              # 与 HiSTREAM、SCOPE 的对比脚本
├── cases/                  # 输入案例及运行结果（结果可重新生成）
├── gui/                    # TiRT 统一页面、时间序列和角度页面
│   ├── server.py           # TiRT GUI 服务（默认 8765）
│   ├── tirteb_server.py    # 能量平衡 API 兼容模块
├── desktop/                # Electron 桌面窗口和服务生命周期管理
├── package.json            # Electron 开发/打包脚本
├── electron-builder.yml    # Windows 安装包配置
├── plot/                   # 输出 CSV 的方向性绘图工具
├── analysis/               # 结果分析脚本和分析产物
├── scripts/                # HiSTREAM 验证复算、打包和截图脚本
├── docs/                   # 操作说明、理论说明、测试文档与 HiSTREAM 验证包
├── .vscode/                # 可选的 VS Code 运行和调试配置
└── *.{bat,command}         # 启动和打包脚本
```

### VS Code 使用

用 VS Code 打开项目根目录：

```powershell
code C:\work\tirt
```

选择 Python 解释器：

```text
Ctrl+Shift+P
Python: Select Interpreter
```

选择已安装依赖的 Python 或 conda 环境。

项目已提供 `.vscode/settings.json` 和 `.vscode/launch.json`，会自动把项目根目录和上一级目录加入 `PYTHONPATH`。这样从子目录运行脚本时也能找到：

```python
from base.xxx import *
```

### 常见问题

#### No module named 'base'

原因是 Python 搜索路径没有包含项目根目录。请从 `C:\work\tirt` 运行，或设置：

```powershell
$env:PYTHONPATH="C:\work\tirt;C:\work"
```

#### No module named 'scipy'

当前环境没有安装依赖：

```powershell
pip install -r requirements.txt
```

#### planck() missing argument

新版 `planck` 已兼容两种写法：

```python
planck(temperature)
planck(wavelength, temperature)
```

建议新代码优先显式传入波长：

```python
planck(10.5, temperature)
```

## 二、TiRT 角度验证（与 HiSTREAM 三维路径追踪对比）

光学 670/860 nm 的 10 场景、1,204 条有效配对见[光学验证报告](docs/TIRT_光学验证报告.md)，整体 RMSE 分别为 0.005374、0.013891 BRF，包含逐场景曲线及 LAI 多光谱、多角度示例。

热红外原主页的 11 场景验证见[热红外验证报告](docs/TiRT_validation.md)。后续 2026-10-21 模型复算另见[版本验证记录](docs/validation/current_model_20261021/README.md)，与原主页报告分别保留。

## 三、界面介绍

### TiRT 桌面窗口

能量平衡、时间序列和多角度计算已经属于 TiRT 的功能模块，用户只启动一个 TiRT 桌面窗口：

```bash
npm install
npm run desktop
```

TiRT 窗口使用 Electron 提供，不会显示浏览器标签栏、地址栏或浏览器菜单；页面中的 `Energy balance` 进入时间序列能量平衡，`Directional model` 返回方向性计算。两者共用同一个服务和项目数据，关闭窗口会同时关闭后台服务。

底层批处理可直接调用 `optical_run.py`、`thermal_run.py` 和 `tirteb_run.py`，不需要启动 GUI 服务；`run.py` 保留兼容。

### 本地页面（Three.js 界面）

`gui/` 提供 TiRT 页面，用于设置输入、预览地表场景并查看方向性和时间序列结果：

```bash
python gui/server.py --port 8765
```

打开 `http://127.0.0.1:8765/`。界面中的 Run simulation 会调用根目录现有的 `run()`，结果区域支持极坐标图、太阳主平面图和垂直太阳主平面图；也可以保存当前表单设置的 `input.csv`。能量平衡页面为 `http://127.0.0.1:8765/time.html`，操作步骤见 [能量平衡操作说明](docs/TIRTEB_user_guide.md)。

### 本地一键启动

项目包含 `environment.yml`、macOS 启动脚本 `start_tirt.command` 和 Windows 启动脚本 `start_tirt.bat`。将整个 `tirt` 文件夹复制到另一台已安装 Miniconda 的电脑后，首次运行启动脚本会自动创建 `tirt` 环境并安装依赖，之后每次运行即可启动本地界面。脚本只绑定 `127.0.0.1`，不会对外提供服务。

如果 macOS 不允许直接双击运行，可在项目目录执行一次：

```bash
chmod +x start_tirt.command
./start_tirt.command
```

Windows 用户请安装 Miniconda，然后双击 `start_tirt.bat`。脚本首次运行会根据 `environment.yml` 创建 `tirt` 环境并安装 Python 依赖，同时安装 Electron 桌面依赖，之后启动无浏览器栏的 TiRT 窗口。也可以在 Anaconda Prompt 中执行：

```bat
cd C:\work\tirt
start_tirt.bat
```

### 打包发布

当前项目提供基于 Python/Conda 的开发运行版本和 Electron 桌面发布版本。Windows 下运行 `python scripts/build_installer.py` 依次生成 PDF 文档、代码压缩包、后台服务（PyInstaller）和最终安装包 `release/v1.2.0/TiRT-Setup-1.2.0-x64.exe`；安装后的 `TiRT` 同时包含方向性和能量平衡时间序列窗口。macOS 上不能可靠构建 Windows `.exe`。

## 参考文献、联系方式与许可证

### 参考文献

#### Physical

1. Zunjian Bian, Shengbiao Wu, Jean-Louis Roujean, Biao Cao, Hua Li, Gaofei Yin, Yongming Du, Qing Xiao, Qinhuo Liu. A TIR forest reflectance and transmittance (FRT) model for directional temperatures with structural and thermal stratification. Remote Sensing of Environment.
2. Zunjian Bian, Biao Cao, Hua Li, Yongming Du, Wenjie Fan, Qing Xiao, Qinhuo Liu. The Effects of Tree Trunks on the Directional Emissivity and Brightness Temperatures of a Leaf-Off Forest Using a Geometric Optical Model. IEEE Transactions on Geoscience and Remote Sensing, 2020.

#### Semi-Physical

3. Zunjian Bian, Biao Cao, Hua Li, Yongming Du, Jean-Pierre Lagouarde, Qing Xiao, Qinhuo Liu. An Analytical Four-component Directional Brightness Temperature Model for Crop and Forest Canopies. Remote Sensing of Environment, 2018.
4. Zunjian Bian, Qing Xiao, Biao Cao, Yongming Du, Hua Li, Heshun Wang, Qiang Liu, Qinhuo Liu. Retrieval of Leaf, Sunlit Soil, and Shaded Soil Component Temperatures Using Airborne Thermal Infrared Multiangle Observations. IEEE Transactions on Geoscience and Remote Sensing, 2016.

#### Semi-Empirical

5. Zunjian Bian, J. L. Roujean, J. P. Lagouarde, Biao Cao, Hua Li, Yongming Du, Qiang Liu, Qing Xiao, Qinhuo Liu. A semi-empirical approach for modeling the vegetation thermal infrared directional anisotropy of canopies based on using vegetation indices. ISPRS Journal of Photogrammetry and Remote Sensing, 2020.
6. Zunjian Bian, Jean-Louis Roujean, Biao Cao, Yongming Du, Hua Li, Philippe Gamet, Junyong Fang, Qing Xiao, Qinhuo Liu. Modeling the directional anisotropy of fine-scale TIR emissions over tree and crop canopies based on UAV measurements. Remote Sensing of Environment, 2021.

### 联系方式

Email: bianzj@aircas.ac.cn

### 许可证 / License

本项目采用 [PolyForm Noncommercial License 1.0.0](LICENSE)。允许依照许可证进行非商业使用、修改和分发；任何商业使用均须事先取得许可方的单独书面授权。

This project is licensed under the [PolyForm Noncommercial License 1.0.0](LICENSE). Noncommercial use, modification, and distribution are permitted under its terms. Commercial use requires a separate prior written license from the licensor.
