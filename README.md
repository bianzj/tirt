# TiRT

Thermal infrared radiative transfer model for vegetation, terrain, urban surfaces, and coupled scenes.

TiRT 用于典型复杂地表热红外方向性辐射传输建模，可计算不同观测角度下的方向亮温、辐射亮度和有效发射率。

## 统一运行入口

项目根目录提供了与 `tirteb` 类似的单点运行方式。默认读取根目录的 `input.csv`，结果写入 `cases/<case_name>/output.csv`：

```bash
python run.py
```

也可以在 Python 中调用外层函数，或通过命令行切换模型：

```python
from run import run

output_path = run("input.csv")
print(output_path)
```

```bash
python run.py --input input.csv --output cases/row_demo/output.csv --model row
```

统一入口支持植被与地表的组合模型：`plane`、`terrain`、`urban` 地表可分别组合 `bare`、`hom`、`row`、`crown` 植被。参数集中在 `input.csv`，输出按“观测方向 × 波段”记录亮温、辐亮度和各组分发射率。

## input.csv 结构

输入分组与 `tirteb` 保持一致：

- `[surface]`：`surface_model` 选择平面、地形或建筑地表，`vegetation_model` 选择裸地、均质植被、垄行或森林。
- `[vegetation_structure]`、`[row_structure]`、`[crown_structure]`：分别配置 LAI、热点、垄行和树冠结构。
- `[terrain]`、`[urban]`：配置地形层和建筑形状，可用分号提供多个高度、半径、密度或建筑参数。
- `[spectral]`、`[urban_spectrum]`：配置热红外高光谱波段，以及叶片、土壤、屋顶、墙壁和街道发射率；每组常数发射率的数量必须与 `wavelengths` 数量一致。
- `[thermal]`：配置土壤、叶片、地形、屋顶、墙壁和街道的日照/阴影温度。
- `[geometry]`：可直接填写角度列表，也可使用 `geometry_source=file` 读取文件。

这里的光谱是热红外高光谱，不是可见光/近红外反射率。光谱文件支持：叶片文件为 `wavelength leaf_emissivity`，土壤文件为 `wavelength soil_emissivity`，建筑文件为 `wavelength roof_emissivity wall_emissivity street_emissivity`；也兼容 `tirteb` 的 7 列建筑文件，读取其中最后三列发射率。文件中的波段会插值到 `wavelengths`。

几何文件每行支持 `vza vaa`，也支持 `vza vaa sza saa`；角度单位为度，注释行以 `#` 开头。

## 适用场景

- 均质植被
- 垄行作物
- 离散树冠
- 坡面植被
- 复杂地形
- 城市建筑
- 城市建筑与街道植被耦合场景

## 方向性绘图

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

## 目录结构

```text
tirt/
├── base/                # 全部现行模型、工具和示例入口
│   ├── hom.py           # 均质植被模型
│   ├── row.py           # 垄行作物模型
│   ├── crown.py         # 离散树冠模型
│   ├── *_voxel.py       # 体素模型
│   ├── terrain*.py      # 地形模型
│   ├── urban*.py        # 城市模型
│   └── run_*.py         # 示例入口
├── plot/                # 极坐标图和太阳主平面图
└── old/                 # 历史代码
```

## 环境安装

建议使用 Python 3.10+。

```powershell
cd C:\work\tirt
py -m venv .venv
.\.venv\Scripts\activate
pip install -r requirements.txt
```

如果 `GDAL` 安装失败，可先安装核心依赖运行主要模型：

```powershell
pip install numpy scipy matplotlib pandas scikit-learn seaborn xlrd opencv-python netCDF4 pysolar
```

## VS Code 使用

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

如果仍出现 `ModuleNotFoundError: No module named 'base'`，请重启 VS Code 终端，或在终端手动设置：

```powershell
$env:PYTHONPATH="C:\work\tirt;C:\work"
```

## 运行示例

```powershell
python -m base.run_tirt_veg
python -m base.run_tirt_vegvoxel
python -m base.run_tirt_slopeveg
python -m base.run_tirt_terrain
python -m base.run_tirt_terrain_veg
python -m base.run_tirt_urban_bt
python -m base.run_tirt_urbanveg_bt_polar
```

## 基本用法

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

## 核心模型

- `Hom`: 均质植被解析模型
- `Row`: 垄行作物解析模型
- `Crown`: 离散树冠解析模型
- `Hom_Voxel`: 均质植被体素模型
- `Row_Voxel`: 垄行作物体素模型
- `Crown_Voxel`: 离散树冠体素模型
- `Slope_Veg`: 坡面植被模型
- `Terrain`: 地形模型
- `Terrain_Veg`: 地形-植被耦合模型
- `Urban`: 城市建筑模型
- `Urban_Veg`: 城市-植被耦合模型

## 常见问题

### No module named 'base'

原因是 Python 搜索路径没有包含项目根目录。请从 `C:\work\tirt` 运行，或设置：

```powershell
$env:PYTHONPATH="C:\work\tirt;C:\work"
```

### No module named 'scipy'

当前环境没有安装依赖：

```powershell
pip install -r requirements.txt
```

### planck() missing argument

新版 `planck` 已兼容两种写法：

```python
planck(temperature)
planck(wavelength, temperature)
```

建议新代码优先显式传入波长：

```python
planck(10.5, temperature)
```

## 参考文献

### Physical

1. Zunjian Bian, Shengbiao Wu, Jean-Louis Roujean, Biao Cao, Hua Li, Gaofei Yin, Yongming Du, Qing Xiao, Qinhuo Liu. A TIR forest reflectance and transmittance (FRT) model for directional temperatures with structural and thermal stratification. Remote Sensing of Environment.
2. Zunjian Bian, Biao Cao, Hua Li, Yongming Du, Wenjie Fan, Qing Xiao, Qinhuo Liu. The Effects of Tree Trunks on the Directional Emissivity and Brightness Temperatures of a Leaf-Off Forest Using a Geometric Optical Model. IEEE Transactions on Geoscience and Remote Sensing, 2020.

### Semi-Physical

3. Zunjian Bian, Biao Cao, Hua Li, Yongming Du, Jean-Pierre Lagouarde, Qing Xiao, Qinhuo Liu. An Analytical Four-component Directional Brightness Temperature Model for Crop and Forest Canopies. Remote Sensing of Environment, 2018.
4. Zunjian Bian, Qing Xiao, Biao Cao, Yongming Du, Hua Li, Heshun Wang, Qiang Liu, Qinhuo Liu. Retrieval of Leaf, Sunlit Soil, and Shaded Soil Component Temperatures Using Airborne Thermal Infrared Multiangle Observations. IEEE Transactions on Geoscience and Remote Sensing, 2016.

### Semi-Empirical

5. Zunjian Bian, J. L. Roujean, J. P. Lagouarde, Biao Cao, Hua Li, Yongming Du, Qiang Liu, Qing Xiao, Qinhuo Liu. A semi-empirical approach for modeling the vegetation thermal infrared directional anisotropy of canopies based on using vegetation indices. ISPRS Journal of Photogrammetry and Remote Sensing, 2020.
6. Zunjian Bian, Jean-Louis Roujean, Biao Cao, Yongming Du, Hua Li, Philippe Gamet, Junyong Fang, Qing Xiao, Qinhuo Liu. Modeling the directional anisotropy of fine-scale TIR emissions over tree and crop canopies based on UAV measurements. Remote Sensing of Environment, 2021.

## 联系方式

Email: bianzj@aircas.ac.cn
