# TIRTEB 理论与方法说明

TIRTEB 是单点、多时相的地表辐射传输与能量平衡模型。每个气象时间步使用太阳几何、地表结构、材料光谱、空气动力学阻抗和植被生理参数，计算组件温度、能量通量及不同方向的观测量。

## 计算流程

~~~text
input_tirteb.csv / meteo.txt / spectral files / geometry files
  -> local solar geometry and surface structure
  -> shortwave and longwave radiation transfer
  -> aerodynamic and surface resistance
  -> photosynthesis or Ball-Berry stomatal resistance
  -> LE, H, G and component energy balance
  -> state update and directional observations (thermal infrared via TiRT)
  -> output.csv, observations.csv, intermediate.csv
~~~

相邻时间步继承前一步的土壤温度、组件温度、热历史和深度温度剖面。能量平衡迭代达到收敛条件或 max_iter 上限后进入下一时间步。

## 地表和植被

- plane：平面地表。
- slope：单坡，`slope_angle`、`slope_aspect`（坡面朝向，自北顺时针）；能量平衡和观测都在坡面局地角度下计算。
- terrain：复合地形，由 `terrain_height`、`terrain_radius`、`terrain_density`（每平方米圆锥数，总占据率须小于 1）描述的圆锥群。
- urban：屋顶、墙面、街道和空气组成的城市场景，建筑由 `building_density`、`building_length/width/height/azimuth` 描述。
- bare：裸土；hom：均质冠层；row：垄行冠层；crown：树冠/森林冠层。复合地形和城市上的植被一律按离散树冠（crown）计算，与 TiRT 相同。

城市和植被可以组合。城市建筑几何独立决定空气动力学阻抗；街道树冠由 `crown_density`（每平方米街道的树冠数）、`crown_rad_a`、`crown_rad_b`、`crown_height` 和 `lai`（每平方米街道）描述。能量平衡中街道被树冠遮挡的比例取 TiRT 树冠模型在天底方向的叶片可视比例，方向反射率使用对应观测方向的叶片可视比例。

## 参数体系

TIRTEB 与 TiRT 使用同一套参数名称和含义（`input_tirteb.csv` 与 `input.csv`）。两个模板共有的参数由同一组读取函数（`base/inputio.py`）解析：

| 类别 | 参数 |
| --- | --- |
| 地表与植被 | `surface_model`、`vegetation_model`、`lai`、`G`、`hspot`（叶片尺寸/冠层高度） |
| 坡面与地形 | `slope_angle`、`slope_aspect`、`terrain_height`、`terrain_radius`、`terrain_density` |
| 垄行与树冠 | `row_width`、`row_blank`、`row_height`、`row_azimuth`、`crown_rad_a`、`crown_rad_b`、`crown_height`、`crown_density` |
| 城市 | `building_density`、`building_length`、`building_width`、`building_height`、`building_azimuth` |
| 热红外发射率 | `emissivity_leaf`、`emissivity_soil`、`emissivity_terrain`、`emissivity_roof`、`emissivity_wall`、`emissivity_street`，以及 `*_spectrum_source/file` 光谱文件 |
| Hapke 土壤 | `hapke_enable`、`hapke_K` |
| 观测几何与波段 | `geometry_mode`、`vza`、`vaa`、`geometry_file`、`geometry_pair_mode`、`principal_vza*`、`spectral_source`、`wavelengths`、`spectrals_file` |

只在 TIRTEB 中出现的是能量平衡参数：气象与时间、`canopy_height`（均质冠层的空气动力学高度；垄行用 `row_height`，树冠用 `crown_height + 2·crown_rad_b`）、`leaf_width`、空气动力学、生理、土壤和城市热参数，以及短波光学参数（`leaf_optics_*`、`soil_optics_*`、PROSPECT、BSM、`roof/wall/street_optics_file`、`reflectance_*`、`absorption_*`）。只在 TiRT 中出现的是组分温度（`*_sunlit_c`、`*_shaded_c`），在 TIRTEB 中它们由能量平衡求出。

## 光谱与观测

可见光和近红外波段使用叶片 PROSPECT 与土壤 BSM 光谱，也可以用 `leaf_optics_file`、`soil_optics_file` 读取实测光谱。热红外直接输入发射率，每个热红外波段（≥ 2.5 µm）一个值或所有波段共用一个值。能量平衡使用各热红外波段发射率的平均值，叶片和土壤的热红外反射率为 `1 - 发射率`（叶片热红外透过率取 0）。`soil_tir_water_effect=1` 时，土壤发射率随含水量从 `emissivity_soil` 线性增加到 `soil_tir_emissivity_wet`。

热红外方向辐亮度和亮温直接调用 TiRT 的辐射传输（`run._surface_result`），输入为能量平衡求得的组分温度：平面用 Hom/Row/Crown 模型，单坡用坡面局地几何（树冠用 Slope_Veg），复合地形用 Terrain/Terrain_Veg，城市用 Urban/Urban_Veg（含树冠与背景之间的散射）。太阳在地平线以下或位于坡面背面时，所有组分按阴影温度计算。因此同一份输入在给定相同组分温度时，TIRTEB 与 TiRT 的热红外结果完全一致（`tests/test_tirteb_tirt_consistency.py`）。observations.csv 对小于 2.5 µm 的波段保存反射率，对热红外波段默认保存亮温 K；坡面背面不可见的方向记为空值。

观测方向由 `geometry_mode` 选择：0 为 `vza/vaa`（`geometry_pair_mode=cross` 为笛卡尔积，`zip` 为一一对应），1 为太阳主平面和垂直太阳主平面，2 为 `geometry_file`，3 为主平面加附加方向。模式 1 和 3 的主平面按每个时间步的实际太阳方位生成。观测时间支持全部时间、精确 DOY 和日内小时分数筛选。

## 能量平衡

对土壤、叶片和城市组件，模型计算净辐射、潜热 LE、感热 H 与热储存 G。能量残差可在 intermediate.csv 中检查：

~~~text
residual = net_radiation - LE - H - G
~~~

土壤热量方法支持 direct、surface 和 depth；气孔阻抗支持 Farquhar 和 Ball-Berry；温度修正支持 direct 和 robust。城市建筑组件使用 `urban_rss`，街道冠层叶片使用 Farquhar 或 Ball-Berry 表面阻抗；城市建筑的潜热为 0，植被潜热按覆盖率聚合。

### 组分通量与场景通量

与 SCOPE ebal（Van der Tol et al. 2009）相同，叶片通量按单位叶面积、土壤通量按单位土壤面积计算，再换算到单位地面：

~~~text
场景通量 = Ps·F土壤光照 + (1 − Ps)·F土壤阴影 + LAIsun·F叶片光照 + LAIshade·F叶片阴影
LAIsun = (1 − Ps) · cos(sza) / G，LAIshade = LAI − LAIsun，G = 0.5
~~~

其中 Ps 是从太阳方向看到的土壤比例（来自 TiRT 冠层模型，垄行和树冠同样适用），夜间 Ps = 0、LAIsun = 0。output.csv 的 Rn、LE、H、G 均为单位地面通量。

组分辐射输入：
- 光照叶片的直射：水平直射辐照度 × G/cos(sza) × 叶片吸收率（每单位叶面积）。
- 叶片的散射光：SAIL 四流求冠层总吸收（入射 − 冠层顶反射 − 土壤吸收），扣除未散射直射后平均到每单位叶面积；土壤得到冠层底的直射和漫射。短波在冠层、土壤和反射之间守恒。
- 长波：冠层分层（每层 ΔLAI ≤ 0.1）求上、下行长波，每层光照叶比例为 exp(−G·x/cos(sza))；叶片两面吸收 ε(E↓ + E↑)，发射 2εσT⁴，土壤吸收 εs·E↓(底)。不再让每片叶片都直接面对天空。
- 叶片和土壤的空气动力学阻抗按 SCOPE 取 (LAI + 1)·(raa + rawc) 和 (LAI + 1)·(raa + raws)。
- 光合与气孔：光照叶和阴影叶分别用各自的叶温和光照（阴影叶只有漫射 PAR）；叶面水汽压按 SCOPE heatfluxes 取 `ea + (ei − ea)·ra/(ra + rs)`，其中 ea 为气象数据的实际水汽压。

`benchmark/scope/` 用 SCOPE 1.73 和 SCOPE 2 自带的验证算例对比上述方案。

### 数值求解与物理方案

- 组件温度用牛顿迭代求解，导数为完整雅可比：`dR/dT = f·4εσT³ + ρcp/ra + ρλ(0.622/p)s/(ra+rs) + dG/dT`，叶片双面发射取 f=2；`dG/dT` 随热通量方法变化（surface 为 `2Γ/(√π√dt)`，depth 为 `λ/dz`，建筑围护结构为 `λ·layers/thickness`）。松弛系数 0.9，robust 模式在残差变号时减半松弛。
- 长波发射和吸收都使用组件发射率（基尔霍夫一致）。
- Monin–Obukhov 稳定度：Obukhov 长度由场景感热 `L = −ρcp u*³ Ta /(κ g H)` 更新（|L| ≥ 1 m），不稳定层结用 Paulson (1970)/Businger–Dyer，稳定层结用 Webb (1970) `ψ = −5ζ`（ζ ≤ 1）。L 作为外层迭代，在每个 L 下能量平衡闭合后再用区间割线法更新。
- 低风速保护（Zeng et al. 1998/CLM）：`U = sqrt(max(u, 1 m/s)² + w*²)`，`w*` 为对流速度尺度（混合层高度 1000 m）。
- 裸土空气动力学阻抗：`ra = [ln(z/z0m) − ψm][ln(z/z0h) − ψh]/(κ²U)`，`z0h = z0m/10`。
- 土壤表面阻抗默认随含水量变化（Sellers et al. 1992）：`rss = exp(8.206 − 4.255·W/Wsat)`；`soil_rss_method=constant` 时使用输入的 `soil_rss`。
- 地表热通量 surface 方法的卷积记忆为 48 小时；depth 方法用第一层土壤节点计算梯度，土壤剖面采用隐式（向后欧拉）导热更新，扩散率为 `λ/(ρ·cs)`。
- 首个时间步长由 `dtime`（秒，默认 1800）给出，之后由气象时间差决定。

城市空气动力学阻抗采用统一的动量/热量 bulk-transfer 形式。屋顶以气象参考风速和城市冠层上方的位移高度、粗糙度及稳定度计算；墙面和街道使用建筑冠层内的局地街谷风速及由建筑平面尺度确定的街谷粗糙度。街谷阻抗设置相对屋顶的物理下限，以体现建筑遮蔽造成的弱交换。街道内叶片的 `raa_leaf` 在街谷阻抗基础上叠加由 LAI、叶片宽度和局地叶片风速决定的叶片边界层阻抗。

因此城市场景中不再令 `raa_leaf` 直接等于墙面阻抗，且 `raa_soil` 与街道表面使用同一街谷交换阻抗。典型排序为 `raa_roof < raa_street = raa_wall < raa_leaf`；实际数值仍随建筑几何、风速和稳定度变化。

## 输出

- output.csv：每个气象时间步的空气、土壤/叶片温度、Rn/LE/H/G（单位地面）、垂直亮温、红光和近红外反射率。
- observations.csv：选定时间、方向和波段的方向性结果。
- intermediate.csv：启用 output_intermediate=yes 后输出的迭代、阻抗、组件比例、`lai_sunlit/lai_shaded/soil_sunlit_fraction` 通量权重和能量诊断。

当前实现面向单点时间序列，不是空间栅格模型。论文或方法比较应同时保存有效的 `input_tirteb.csv`、`output.csv`、`observations.csv` 和 `intermediate.csv`。
