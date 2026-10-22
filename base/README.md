# `base/` 核心算法

`base/` 是 TiRT 和 TIRTEB 共用的 Python 核心包。根目录入口负责参数和输出，模型计算尽量放在这里，避免 GUI 或脚本复制算法。

## 模块分层

- `hom.py`、`row.py`、`crown.py`：均质植被、垄行和离散树冠模型。
- `slope_veg.py`、`terrain*.py`、`urban*.py`：坡面、复杂地形和城市场景（地形植被由 `terrain_forest.py` 计算）。
- `gap.py`、`path.py`、`hotspot.py`、`proportion.py`、`scatter.py`：几何路径、间隙、热点、方向比例和散射计算。
- `radiance.py`、`radiation.py`、`rt_transfer.py`：辐射、亮温和辐射传输公共逻辑。
- `aerodynamic.py`、`biochemical.py`、`evapo.py`、`eb.py`：TIRTEB 的阻抗、蒸散和能量平衡。
- `io.py`、`inputio.py`、`modelio.py`、`structs.py`：输入、数据文件和模型结构解析。

外部代码应从根目录启动，并通过 `from base...` 导入。不要在 `base/` 中写入运行结果、临时文件或 GUI 代码。

## 参数化标准与 TIRTEB 衔接

`crown_parametric.py`、`row_parametric.py` 为默认的树冠、垄行闭式四组分比例；`thermal_crown.py`、`thermal_row.py` 保留此前的 `clumping_dual`、`conditioned_corridor` 可选分支；`terrain_forest.py`、`urban_veg.py` 为地形与城市植被；`tirteb_tirt.py` 把 TIRTEB 的参数和组分温度映射到 `run._surface_result`；`single_point.py` 是单场景能量平衡入口。公式见[理论说明](../docs/TiRT_theory_and_methods.md)，与 HiSTREAM 的对比见[测试文档](../docs/TiRT_validation.md)。
