# `base/` 核心算法

`base/` 是 TiRT 和 TIRTEB 共用的 Python 核心包。根目录入口负责参数和输出，模型计算尽量放在这里，避免 GUI 或脚本复制算法。

## 模块分层

- `hom.py`、`row.py`、`crown.py`：均质植被、垄行和离散树冠模型。
- `*_voxel.py`：对应的体素模型。
- `slope_veg.py`、`terrain*.py`、`urban*.py`：坡面、复杂地形和城市场景。
- `gap.py`、`path.py`、`hotspot.py`、`proportion.py`、`scatter.py`：几何路径、间隙、热点、方向比例和散射计算。
- `radiance.py`、`radiation.py`、`rt.py`、`rt_transfer.py`：辐射、亮温和辐射传输公共逻辑。
- `aerodynamic.py`、`biochemical.py`、`evapo.py`、`eb.py`：TIRTEB 的阻抗、蒸散和能量平衡。
- `io.py`、`inputio.py`、`fileio.py`、`modelio.py`、`structs.py`：输入、数据文件和模型结构解析。
- `run_tirt_*.py`：历史低层示例和回归入口，供复现实验使用；新代码优先调用根目录的 `run.py` 或 `tirteb_run.py`。

外部代码应从根目录启动，并通过 `from base...` 导入。不要在 `base/` 中写入运行结果、临时文件或 GUI 代码。
