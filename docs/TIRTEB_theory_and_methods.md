> 当前标准版的完整公式见 [TiRT/TIRTEB 理论与方法说明](TiRT_TIRTEB_理论与方法说明.pdf)（源文件 `docs/report/TIRT_TIRTEB_理论与方法.md`）。本文件包含 legacy/upgrade 等历史选项，仅作参考。

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

相邻时间步继承组件温度、热历史、深度温度剖面和上一步稳定度。每个时间步先在给定 Obukhov 长度下求组件能量平衡，再用场景感热更新稳定度；达到联合闭合条件或外层迭代上限后输出结果。达到迭代上限不等于收敛，应检查 `closed` 和残差。

## 地表和植被

- plane：平面地表。
- slope：单坡，`terrain_slope`、`terrain_aspect`（坡面朝向，自北顺时针）。输入文件也接受 TiRT 的 `slope_angle`、`slope_aspect` 作为同义键；两种名称同时给出且数值不同时报错。输入气象短波按水平面定义，投影后在局地坡面计算，通量按单位坡面面积输出。
- terrain：复合地形，由 `terrain_height`、`terrain_radius`、`terrain_density`（每平方米圆锥数，总占据率须小于 1）描述的圆锥群。
- urban：屋顶、墙面、街道和空气组成的城市场景，建筑由 `building_density`、`building_length/width/height/azimuth` 描述。
- bare：裸土；hom：均质冠层；row：垄行冠层；crown：树冠/森林冠层。复合地形和城市上的植被由所选冠层模型（hom、row 或 crown）直接计算，与 TiRT 相同。

城市和植被可以组合。城市建筑几何决定屋顶和街谷的空气动力学阻抗，街道植被由所选冠层模型描述。能量预算保留全部街道土壤面积，并另外计入叶面积；天底叶片可视比例只属于方向观测，不能用来扣减街道能量面积。复合地形和城市总通量按单位水平地面面积输出。

## 参数体系

TIRTEB 与 TiRT 尽量使用同名参数；TIRTEB 的字段由 `base/io.py` 读取，TiRT 由 `base/inputio.py` 读取。下表为两者含义相同的参数（TIRTEB 读取时括号中的名称为同义键）：

| 类别 | 参数 |
| --- | --- |
| 地表与植被 | `surface_model`、`vegetation_model`、`lai`、`hspot`（热点参数，无量纲） |
| 坡面与地形 | `terrain_slope`（=`slope_angle`）、`terrain_aspect`（=`slope_aspect`）、`terrain_height`、`terrain_radius`、`terrain_density` |
| 垄行与树冠 | `row_width`、`row_blank`、`row_height`、`row_azimuth`、`crown_rad_a`（水平半径 m）、`crown_rad_b`（垂直半轴 m）、`crown_height`（冠底离地高度 m）、`crown_density`（每平方米树冠数） |
| 城市 | `building_density`、`building_length`、`building_width`、`building_height`、`building_azimuth` |
| 波段参数表 | `wavelengths`、`leaf_reflectance`、`leaf_transmittance`、`soil_reflectance`、`street_reflectance`、`wall_reflectance`、`roof_reflectance`（与 TiRT 同一格式；TIRTEB 允许 `-` 表示取模型/实测值） |
| Hapke 土壤 | `hapke_enable`、`hapke_K` |
| 观测方向与波段 | `vza`、`vaa`、`spectral_source`、`spectrals_file` |

树冠为 Poisson 分布的椭球体，冠层开放度 `exp(−crown_density·π·crown_rad_a²)`，允许树冠相互重叠，因此 `crown_density·π·crown_rad_a²` 可以大于 1。旧的 `crown_radius`、`crown_width` 已停用，出现时报错（旧版 `crown_height` 表示垂直半轴，与现在的冠底高度含义不同，静默兼容会算错）。

观测方向与 TiRT 不同：TIRTEB 用 `direction_source=direct` 时取 `vza` 与 `vaa` 列表的笛卡尔积，`direction_source=file` 时读取 `directions_file`（每行 `vza vaa`）。TiRT 的 `geometry_mode`、`geometry_file`、`geometry_pair_mode`、`principal_vza*` 在 TIRTEB 中不支持，出现时报错。不认识或不起作用的参数名读取时给出警告并忽略。

TIRTEB 还包含气象与时间、`leaf_width`、空气动力学、生理和热参数，以及驱动短波能量预算的叶片、土壤和城市材料光谱。TiRT 中外部给定的组分温度，在 TIRTEB 中由能量平衡求解；土壤与城市表面的初始温度取模拟首日气温的日平均，不由输入给定。参数读取与映射以 `base/io.py`、`base/structs.py` 和当前输入模板为准。

## 光谱与观测

观测波段和各组分光学参数用与 TiRT 相同的 7 行波段参数表给出（波段、叶片反射率、叶片透过率、土壤、街道、墙壁、屋顶反射率，分号分隔，每个波段一个值）。TIRTEB 中任一格可以填 `-`：

- 短波波段（< 2.5 µm）：填了数值就以表格为准，用于该波段的观测反射率；填 `-` 时取叶片 PROSPECT、土壤 BSM（或 `vegetation_spectrum_file`、`soil_spectrum_file` 实测光谱）以及城市光谱在该波段的值。能量平衡的全波段短波吸收始终使用 PROSPECT/BSM 或实测光谱，不受表格短波值影响。
- 热红外波段（≥ 2.5 µm）：按基尔霍夫定律换算发射率，叶片 `ε = 1 − ρ − τ`，土壤、街道、墙壁、屋顶 `ε = 1 − ρ`；填 `-` 时城市组分取实测城市光谱中的发射率（`urban_spectrum_source=measured`），否则取与 TiRT 相同的默认值（叶片 ρ 0.015、τ 0、土壤 0.045、街道 0.045、墙壁 0.07、屋顶 0.05）。观测表仅决定对应观测波段的光学性质；能量平衡使用独立的宽带灰体参数，不能由观测波段算术平均决定。`soil_tir_water_effect=1` 时，表中土壤值作为干土发射率，随含水量线性增加到 `soil_tir_emissivity_wet`。

宽带长波默认采用固定 10.5 µm 的灰体参考值，与观测波段数量、间距和取值无关。叶片默认 ε=0.985、τ=0，土壤默认 ε=0.955，屋顶/墙面/街道默认 ε=0.95/0.93/0.955；实测城市发射率光谱可改变其固定 10.5 µm 参考值，启用土壤含水量发射率关系也会改变土壤灰体参考值。这是灰体近似，不是完整热谱的 Planck 加权和半球积分。

以下可选标量直接指定能量平衡参数，不覆盖观测表：`broadband_leaf_emissivity`、`broadband_leaf_transmittance`、`broadband_soil_emissivity`、`broadband_urban_emissivity_roof`、`broadband_urban_emissivity_wall`、`broadband_urban_emissivity_street`。各值必须有限且在 [0,1] 内，叶片发射率与透过率之和不得大于 1。显式宽带值按输入使用。

`spectral_source=file` 时从 `spectrals_file`（默认 [observation_spectrals_tirteb.txt](../data/spectrals/observation_spectrals_tirteb.txt)）读取同样 7 行的表，每行前面可以带文字标签。

土壤短波反射率由 BSM 计算；土壤含水量统一用体积分数 m³ m⁻³，BSM 的容量参数也取 0.25，并直接在该尺度上计算 `μ = (θ − 0.05)/0.25`，不再换算成百分数。为兼容旧输入，`soil_moisture`、`satwater` 和 TIR 干/湿参考含水量中大于 1 的值按百分数除以 100，模型内部只保留 m³ m⁻³；`bsm_latitude`、`bsm_longitude` 是 BSM 的谱形角度参数（度），与土壤含水量无关。水膜反射率取 `1 − tav(40°, nw)`（`calctav` 返回的是透过率），土壤水膜下界反射率按 SCOPE 的 `rbac` 公式计算。

短波默认按直射 0.8、漫射 0.2 划分（`shortwave_diffuse_fraction=0.2`）；填 `-1` 时按 Erbs 从总辐射和太阳高度角反算，填 0–1 时使用给定值。TiRT 方向模型的短波反射率使用同一默认值。

能量平衡默认还启用 `shortwave_enforce_solar_limit=1`：直射法向辐照度不得超过晴空直射辐照度。晴空 DNI 采用 Meinel 模型 `DNI_clear = 1367·E0·0.7^(AM^0.678)`（Meinel & Meinel, 1976），E0 = `1+0.033·cos(2π·DOY/365)`，相对大气质量 AM 采用 Kasten & Young (1989) 公式 `AM = 1/[cos(sza) + 0.50572·(96.07995−sza)^−1.6364]`（sza 为度）。若给定比例下的直射超过这一上限，保留原气象总短波，将超出的直射功率转为漫射；地平线以下的正短波全部作为漫射。日地距离修正采用 [FAO 的年内近似](https://www.fao.org/4/x0490e/x0490e07.htm)。该处理是物理可行性约束，不是实际直散射观测反演；Meinel 值是典型晴空值，高海拔极洁净大气下真实 DNI 可能略高于此上限。中间文件保留请求比例、实际比例、`shortwave_partition_limited`、直射法向功率和上限（`solar_direct_normal_limit_W_m2`）；控制试验可显式设为 0 关闭限制。单纯给定几何的静态 BRF 没有绝对入射功率，不进行这一约束。小于 20 W m⁻² 的正入射也继续参与预算，不再直接丢弃。

所有气象值一律按时间戳时刻的瞬时值使用，原始数据为时段均值时也如此（项目约定），太阳几何在时间戳时刻计算，模型不平移时间。每步记录晴空指数 `clearness_index = Rs/(1361·E0·cos(sza))`（太阳在地平线以下时为 NaN）；`clearness_index_exceeds_one=True` 表示该步气象总短波超过同一时刻的大气层顶水平辐照度。对时段均值数据，这在日出、日落附近是这一约定的预期结果：保留气象总短波，超过晴空直射上限的部分按漫射处理。运行结束时汇总这些时刻及 kt 值，作为提示性告警。

TIRTEB 的方向观测采用 TiRT 标准实现。平面冠层的吸收预算使用 hom3 / row3 / crown2 节点方案，土壤与叶片之间满足互易交换，并对散射求全部阶；城市与地形在此冠层算子外耦合各自的交换边界。方向短波反射率（`green`、`nir` 与 observations.csv 的短波列）调用 `canopy_closure.reflectance`：一次散射使用与热红外相同的四组分可见比例，高阶散射使用同一节点算子。由于方向一次散射含热点相关而节点吸收预算按独立间隙计一次逃逸，方向反射率的半球积分可比吸收预算对应的反照率高约 0.02–0.03（近红外、LAI≈2）；能量平衡使用吸收预算，不受此差异影响。方向算子细节及已知限制见 [TiRT 理论与方法说明](TiRT_theory_and_methods.md) 第 4.6–4.7 节。

方向观测：平面和单坡的短波 BRF/HDRF 直接调用标准版四组分 + 节点辐射度（单坡在坡面局地坐标下计算，垄向同样换算到局地坐标，再以水平面总入射为分母：局地直射反射率乘 `(1−f)·cos i/cos sza`、天空反射率乘 `f·(1+cos β)/2`，权重取自能量平衡的同一坡面投影，与 TiRT 一致；`hapke_enable=1` 时一次项土壤取 Hapke 方向反射）；复合地形调用 TiRT 标准版复合地形光学；城市调用 TiRT 城市光学（街道冠层同为标准版）。`transmittance` 为标准版视场叶片比例 × 叶片透过率。

热红外方向辐亮度和亮温直接调用 TiRT 标准版的辐射传输（`run._surface_result`，经 `base/tirteb_tirt.py` 映射），输入为能量平衡求得的组分温度。参数映射：`terrain_type=1` → TiRT 单坡（`terrain_slope/aspect` → `slope_angle/aspect`），`terrain_type=2` → TiRT 复合地形（锥体占据率 `density·π·radius²` 必须小于 1）；树冠 `crown_rad_b` → TiRT `crown_rad_b`（内部 `height`），`crown_rad_a` → `crown_rad_a`（内部 `width`），`crown_density` → `crown_density`（内部 `stand_density`），`crown_height` → `crown_height`（内部 `crown_base_height`，冠底离地高度）；叶片发射率 1−ρ−τ、土壤 1−ρ_s，城市组分 1−ρ，均取自该热红外波段的波段参数表；地形锥体按土壤温度处理；街道取能量平衡求得的街道温度。太阳在地平线以下或位于坡面背面时，所有组分按阴影温度计算。同一场景、相同组分温度下 TIRTEB 与 `run.py` 的辐亮度一致（plane/slope 的全部植被组合、terrain 的非裸地组合、urban 的全部植被组合；terrain 裸地的锥体在 TIRTEB 中取土壤温度）。observations.csv 对小于 2.5 µm 的波段保存反射率，对热红外波段默认保存亮温 K；`observation_quantity=reflectance` 时热红外保存 1 − 方向发射率（同一 TiRT 模型在等温黑天空下的辐亮度 / B(T)，含全部散射阶）。坡面背面不可见的方向记为空值。方向观测与 TiRT 的一致性可用 `validation/compare_tirt_tirteb.py --full` 复算（相同场景、组分温度和波段光学下逐方向比较，默认容差 `1e-9`；场景用 `input_tirteb.csv` 的参数名经 `read_input/read_defined` 读入，TiRT 参考的几何直接取自这些参数，因此也能检出参数映射错误；算例中太阳均在地平线以上）。

观测方向由 `direction_source` 选择：`direct` 为 `vza` 与 `vaa` 两个列表的笛卡尔积，`file` 为 `directions_file` 中逐行给出的 `vza vaa`；方向在整个时间序列中固定，不随太阳方位生成主平面。观测时间由 `observation_time_mode` 选择：`all` 为全部时间；`time` 为 `observation_times` 中精确的 DOY 时刻；`hour` 为 `observation_hours` 中的当地钟点（0–24 小时，可带小数，例如 `12;13.5`），每天相同钟点都输出，非法值报错。

## 能量平衡与面积基准

各组分能量残差定义为 `R = Rn − LE − H − G`，正 G 表示热量从表面进入材料。裸土和普通植被可蒸发/蒸腾；当前城市分支将屋顶、墙面、街道的 LE 置零，街道冠层叶片计算 LE。因此城市结果不能直接解释为具有湿街道蒸发的水量闭合模型。

### 平面冠层

土壤通量按单位土壤面积计算。叶片采用单位冠层地面的“大叶”方程：输入为 LAI × 单位叶面积吸收，发射为 `2·LAI·εσT⁴`。边界层阻抗包含 `70/LAI·sqrt(leaf_width/u_leaf)`；生理模型先给出单位叶面积的气孔阻抗，再除以 LAI 转换到同一大叶面积基准。`rssunlit_leaf / rsshaded_leaf` 保存转换前的叶面积阻抗。

平面的场景通量为：

~~~text
F = Ps·Fsoil_sun + (1−Ps)·Fsoil_shade
    + (LAIsun/LAI)·Fleaf_sun + (LAIshade/LAI)·Fleaf_shade
LAIsun = 节点截获直射份额 / (0.5/cos(sza))
LAIshade = LAI − LAIsun
~~~

叶面积划分的太阳天顶角计算上限为 85°。夜间直接照明面积为零；面积为零的组分不贡献场景总量。上述叶片权重与大叶方程的 LAI 因子配合，等价于按真实叶面积积分。`soil_sunlit_fraction / soil_shaded_fraction / lai_sunlit / lai_shaded` 保存能量面积权重；复合地形的土壤权重和为真实表面积比，可以大于 1。

平面短波在 2001 个 1 nm 波段上求解。一次直射只给光照组分；直射的后续散射和天空漫射用节点辐射度算子求全部阶，并按节点中光照/阴影面积分配吸收。总叶片吸收、土壤吸收和逃逸反射应与入射相平衡。长波也在同一互易节点几何上求全部阶，节点发射由光照/阴影温度的四次方按面积混合；等温节点处于同温黑体辐射环境时净长波为零。

### 城市

建筑单位水平地面的面积为：

~~~text
Aroof = density·length·width
Awall = density·2·height·(length+width)
Astreet = 1−Aroof
~~~

这些面积不再除以建筑外表面积总和。屋顶、墙面和街道按实际日照份额加权：屋顶在白天接受直射，墙面按建筑方位和各立面朝向的正太阳投影划分；建筑遮阴下街道的日照份额为 `b·Ps`，叶片为 `b·LAIsun/LAI`。其中 b 是统计建筑透光率 `exp(−obstruction)`，其遮挡深度包含建筑密度、高度、长宽、太阳天顶角和相对建筑方位。夜间和正午垂直入射时，垂直墙面的直射照明面积为零。

场景土壤使用完整 Astreet，叶片在同一街道上另外贡献真实叶面积；不会按天底冠层覆盖率排除街道面积。墙面真实面积比可以大于 1，故正确的墙面地面总量也可以超过其单位墙面通量。

`base/urban_energy.py` 使用街道/冠层边界和四个立面的五节点互易视因子。街道/冠层边界的反射率 R 和热源 E 来自同一个局地冠层吸收/发射算子。长波辐射度满足：

~~~text
J = [I − diag(R)·F]⁻¹ · [E + R·Fsky·Lsky]
Hincident = F·J + Fsky·Lsky
~~~

街道土壤与叶片再用同一冠层算子分配这一入射；不会分别套用互不一致的天空/墙面经验混合权重。墙面温度按实际光照面积混合，屋顶对上方天空开放。所有表面和天空等温时，墙面、街道、叶片各自净长波为零。

短波同样耦合这套城市视因子和冠层边界。直接入射功率严格分成屋顶 `Aroof·B`、街道冠层 `Astreet·b·B`、墙面 `Astreet·(1−b)·B`；墙面内部按正太阳投影分配。天空漫射及城市/冠层多次反射随后共同求解。`urban_SW_absorbed + urban_SW_escaped` 等于外部入射；`urban_LW_escaped` 是单位水平地面向外的长波，`urban_LW_net_ground = Lsky − urban_LW_escaped`。

### 单坡与复合地形

气象总短波必须为水平面辐照度，直射 B 与漫射 D 由 `shortwave_diffuse_fraction` 划分。单坡将水平直射投影为 `B·max(cos(i),0)/cos(sza_global)`，背阳坡直射为零；各向同性天空因子为 `(1+cos(slope))/2`。下半球其余视域接受同一材料/冠层参数的周围水平地表反射短波和出射长波，采用当前共享组分温度。单坡的结果按单位局地坡面面积给出，LAI 也按坡面面积定义；外部已订正到坡面的气象短波不能作为水平入射再次输入。

复合地形采用水平裸露地面和 12 个圆锥方位面群体的统计近似。局地太阳角和统计相互遮阴决定首次直射分配，再归一化使穿过单位水平孔口的直射只截获一次；视因子按交换面积对称化，满足互易性和单位水平天空孔口约束。面群体长波通过有效冠层反射/发射边界求全部阶互热。面间短波先由局地谱吸收求有效反射，再以给定漫射谱形回投，属于灰体与谱形结合的近似，不是逐波长保留每次面间反射谱变化的精确解。

复合地形仍只求四个共享温度 Tss/Tsh/Tls/Tlh，没有为每个方位面独立求温度。LAI 按水平地面面积定义，局地面群体 LAI 按坡度余弦换算；所有通量按真实面面积积分后输出到单位水平地面。当前方案不是 DEM 光线追踪或逐坡面温度模型，也没有独立实测地形温度精度验证。检查地形机制应同时看坡面投影、能量面积、辐射闭合及共享温度这一近似。

### 叶片生理与土壤水分

Farquhar 和经验 Ball–Berry 分别使用光照叶 Tls 与阴影叶 Tlh。光照叶 PAR 为直接 PAR 加对应光照叶散射 PAR，阴影叶仅接受对应阴影叶散射 PAR；PAR 单位为 µmol photons m⁻² leaf s⁻¹。温度和光照均相同时两类叶片响应可以相同，否则分别计算。两种方案的气孔导度都不低于 `stomatal_g0`（默认 0.01 mol m⁻² s⁻¹）：Farquhar 方案在净光合 A≤0 时取 g0，与 SCOPE 2 一致，不再按 SCOPE 1.73 的做法把气孔阻抗设为 625000 s m⁻¹（完全关闭）。`stressfactor` 仍是给定的植被胁迫参数，不由新水量状态方程求出。

未显式给定 `soil_rss`/旧别名 `rss` 时，默认使用 Sellers 经验关系：

~~~text
rss = exp(8.206 − 4.255·clip(soil_moisture/satwater,0,1))   [s m⁻¹]
~~~

显式给定 `soil_rss` 时保留常数含义；`soil_rss_method=constant` 强制常数模式，`soil_rss_method=sellers` 强制上述含水量关系。`satwater` 使用实际输入值，`soil_moisture` 为给定体积含水量（m³ m⁻³；大于 1 的百分数输入会除以 100）。模型不在此基础上声称已经闭合降水、渗流、根系吸水和储水变化的水量平衡。经验关系来源：[Sellers et al. (1992)](https://denning.atmos.colostate.edu/biocycle/Documents/SiB/Sellers_1992.pdf)。

## 导热与时间离散

`soil_cs` 是质量比热，单位 J kg⁻¹ K⁻¹；土壤体积热容为 `soil_rhos·soil_cs`，扩散率为 `soil_lambda/(soil_rhos·soil_cs)`。`urban_cs` 本身是体积热容，单位 J m⁻³ K⁻¹，不再乘 `urban_rhos`。模板中的数值必须按这些单位解释。

- direct：G 取给定的净辐射比例（`soil_heat_fraction` 或 `urban_heat_fraction`）。
- surface：用半无限介质表面温度变化的卷积核求 G，热惯量为 `sqrt(lambda·Cvol)`；保存实际经过的 48 小时时间窗及各段真实间隔，包含跨越截断边界的最老线性段。
- depth（土壤及街道）：土柱深 2 m，层厚自表层 0.01 m 起按 1.3 倍几何加粗（共 16 个节点，常量见 `base/heat_conduction.py`），底部为零通量（Neumann）边界。有限体积向后欧拉：节点 i 的热容为 `Cvol·(h_(i−1)+h_i)/2`，底节点为半层，`Cvol·w_i·(Ti−Ti_old)/dt = lambda·(Ti−1−Ti)/h_(i−1) − lambda·(Ti−Ti+1)/h_i`，不截断 Fourier 数。表面节点温度为当前迭代的表面温度，表面热通量包含第一层导热和表层半格储热：`G = lambda·(T0−T1)/h0 + Cvol·(h0/2)·(T0−T0_old)/dt`，因此 `G·dt` 等于整列储热变化。每个气象步只在最终状态提交一次剖面，迭代中不重复推进历史。`soil_depth_step / urban_depth_step` 不再参与土柱网格。正弦表面温度驱动（振幅 10 K）下，第 7–8 天平均 G 约 −0.1 W m⁻²，振幅与解析半无限解相差 <1%；首日平均 G 仍受初始土柱温度（取首日平均气温）与表面平均温度之差的调整过程影响。
- depth（屋顶/墙面）：均匀分层的有限围护结构，内侧为固定室内温度，表面通量取第一层梯度。

首步 `dt=dtime`（秒，默认 1800），之后采用气象时间戳的实际差；时间必须严格递增。城市街道导热与历史推进都使用土壤热参数，屋顶/墙面 depth 使用城市体积热容、有限围护结构、给定厚度/层数和固定室内温度。固定边界、有限层数及 surface 的半无限介质假定各有适用范围。

## 耦合稳定度与数值求解

组件温度首先采用阻尼的分量近似牛顿修正，固定 `WC=0.08`，单次温度变化限制为 ±1.5 K，温度保护范围为 240–380 K。`direct` 使用各分量近似导数；`robust` 使用简化的传热/辐射分母。标量修正达到每批上限仍未闭合时，采用有界信赖域最小二乘后备求解，同时把组分温度和控制阵风速度的场景 H 作为未知量，用数值雅可比求全部组分能量方程及 H 一致性方程。数值扰动期间固定前一时间步的热剖面和历史，不重复推进时间。

傍晚生理响应切换可能改变温度平衡分支，使此前建立的标量 L 区间失效。外层达到上限仍未共同闭合时，再联合求解温度、阵风 H 与 `1/L`，增加方程 `H + rho·cp·Ta·ustar³·(1/L)/(kappa·g) = 0`，最多进行 `4·max_iter` 次信赖域目标求值。后备求解器的 success 标志不替代物理残差检查；最终仍按同一共同闭合判据判断。`stability_joint_evaluations` 记录这一联合后备路径的求值数，`iterations` 包含标量及后备温度求值数。

Monin–Obukhov 外层以 `x=1/L` 为未知量，避免中性点 L 无穷大的奇异表示。每次固定 x 重求组件温度和通量，场景感热 H 使用与输出相同的面积积分（城市含真实墙面积），但稳定度需要单位水平地面通量：单坡的输出 H 按坡面面积计，求 L 及阵风速度时除以 `cos(slope)` 换算到单位水平面积（`stability_flux_area_factor = 1/cos(slope)`，`scene_sensible_heat_W_m2` 即换算后的值）；平面、复合地形和城市本已按水平地面计，因子为 1。L 只用 H，不含 LE 的虚温浮力项。只换算稳定度所用通量，输出通量的面积基准不变。目标为：

~~~text
Ltarget = −rho_air·cp_air·Ta·ustar³ / (kappa·g·H)
xtarget = −kappa·g·H / (rho_air·cp_air·Ta·ustar³)
~~~

空气密度每步由气象计算，`rho_air = p/(Rd·Tv)`，`Tv = Ta/(1−0.378·e/p)`，Rd = 287.05 J kg⁻¹ K⁻¹；p 与水汽压 e 按 hPa 输入并换算为 Pa。H、LE、温度修正分母、阵风速度、L 以及气孔阻抗的 mol m⁻² s⁻¹ ↔ s m⁻¹ 换算（`rho_air/M_air`）均使用这一密度，`air_density_kg_m3` 记录每步取值。

H>0 对应 L<0，H<0 对应 L>0；中性可返回 L=inf，不强行把 L 裁到 ±1 m。只有内层最大组件能量残差达到 1e−4 W m⁻² 后，才将当前点加入稳定度求根的符号区间。区间建立后用有保护的割线更新，并限制在区间内部 10%–90%；区间尚未建立时用 `xnew=(x+xtarget)/2`。内层每批最多 `max_iter` 次（默认 50），外层最多 `max_stability_iter` 次（默认 60）；内层（含后备求解）在当前 L 下仍未闭合时，该点不进入符号区间，但仍以 `xnew=(x+xtarget)/2` 更新 L，避免下一批在同一 L 上重复同一结果。`iterations` 是全部外层批次累计的内层次数，可大于 `max_iter`。

最终返回前重新计算温度对应的通量。`stability_closed` 要求 `abs(1/L−xtarget) ≤ 1e−7 + 1e−6·abs(xtarget) + |c|·1e−4/ustar³`，单位 m⁻¹，其中 `c = −kappa·g/(rho_air·cp·Ta)`，`xtarget = c·H/ustar³`。最后一项是能量容差允许的 H 误差（1e−4 W m⁻²）传递到 1/L 上的量：在 H 接近 0、u* 很小的傍晚和夜间，1/L 不可能比 H 本身算得更准（外层迭代内部使用该容差的一半，使最终复核留有余量）；`closed` 还要求内层闭合、阵风 H 一致性残差不超过 1e−4 W m⁻²，以及最终最大能量残差不超过公共阈值 5 W m⁻²。输出了数值温度或稳定度单项闭合，不代表联合闭合已经成立。未闭合的步照常推进热历史；每步记录 `closed`、`energy_balance_residual_W_m2`（最大组件残差）、`fallback_solver_used`（是否调用了信赖域温度后备或联合求解）及 `stability_limited`。

不稳定层结使用 Businger–Dyer/Paulson 型积分函数。稳定层结先对整个廓线统一限制倒数长度：`xeff=min(1/L,1/zref)`，随后所有高度使用 `psi(z)=−5·max(z,0)·xeff`。这里 zref 为交换的参考高度，不对参考高度和粗糙度高度分别裁剪 z/L，避免强稳定时两端订正相消而恢复中性交换。由 H 求出的实际 L 保持不裁剪。`stability_zeta` 为诊断参考高度除以实际 L；`stability_limited` 由阻抗计算在 `zref/L>1`、即上述截断实际生效时置为 True。该标记说明使用了稳定订正适用范围的上限，不代表强稳定夜间湍流经过独立验证。

有效风速包含对流阵风尺度：

~~~text
wstar = [g·max(H,0)·1000/(rho_air·cp_air·Ta)]^(1/3)
U = sqrt(max(u, min_wind_speed)² + wstar²)
~~~

`u` 为气象文件的原始风速（读入时不截断）；`min_wind_speed` 默认 0.05 m s⁻¹，可在 `input_tirteb.csv` 的 `[atmosphere]` 中修改。

与 SCOPE（`resistances.m`）相同，每一段气动阻抗都不超过 `max_aerodynamic_resistance`（默认 400 s m⁻¹，设为 0 表示不限制）。分段为：冠层上方 rai+rar+rac；叶边界层；冠层内到土壤（含 aero_rbs）；裸土的 bulk-transfer 阻抗；城市的屋顶、街谷 bulk-transfer 及叶边界层阻抗。叶片总阻抗 raa_leaf = 冠层上方段 + 叶边界层段，土壤总阻抗 raa_soil = 冠层上方段 + 冠层内段，所以总阻抗最大可到两倍上限。这是经验限幅：在弱风、稳定层结下防止交换完全停止，不是物理推导。

1000 m 是固定混合层高度近似。裸土使用动量和热粗糙度的 bulk-transfer 阻抗，`z0m=0.01 m`，`z0h=z0m/10`；冠层稳定度修正参与摩擦速度和冠层扩散交换。冠层 d 与 z0 采用 Raupach (1994) 关系（与 SCOPE 相同）：`d = hc·[1−(1−exp(−sqrt(cd1·LAI/2)))/sqrt(cd1·LAI/2)]`，`z0 = (hc−d)·exp(−kappa·G1+psicor)`，`G1 = max(3.3, (cssoil+cr·LAI/2)^−1/2)`，因此稀疏冠层的 z0 随 LAI 减小。为使 LAI→0 时与裸土连续，LAI < 0.5 时以权重 `w = s²(3−2s)`（`s=LAI/0.5`）混合冠层与裸土的 ustar、土壤阻抗 raa_soil、稳定度参考高度，d 按 `w·d`、z0 按对数加权报告（`aerodynamic_displacement / aerodynamic_roughness_length / aerodynamic_canopy_weight`）；LAI ≥ 0.5 时结果与原冠层方案相同。LAI=1e−6、1e−4、0.01 的 raa_soil、ustar 与土壤温度单调趋近裸土结果。城市屋顶采用冠层上方交换，墙面/街道采用局地街谷风速和粗糙度，叶片再加叶边界层阻抗。这些经验传输近似及低风保护不能替代稳定夜间或真实城市的独立通量验证。

## 输出与验证边界

- `output.csv`：组分温度、Rn/LE/H/G、方向观测摘要及联合闭合状态。Rn/LE/H/G 的单位面积基准：平面为地面、单坡为局地坡面、复合地形和城市为水平地面。
- `L`（m）、`ustar`（m s⁻¹）、`stability_closed`、`stability_iterations`、`stability_residual_m_inv`（m⁻¹）和累计 `iterations` 写入 `output.csv`。`energy_balance_residual_W_m2` 是最大组件残差，而非场景平均残差。
- `observations.csv`：选定时间、方向和波段的反射率或亮温。方向观测和面积积分通量属于不同量，不能用天底观测比例替代能量面积。
- `intermediate.csv`：`output_intermediate=yes` 时保存完整状态及 `Rn_total / LE_total / H_total / G_total`；包括各组件 `EB_residual_*`、辐射输入、阻抗、真实面积权重、`stability_zeta / stability_limited / L_target`、热历史、城市辐射闭合等诊断。零面积组分保留有限诊断值，但不贡献场景总量。

`validation/compare_tirt_tirteb.py --full` 检查相同温度和波段光学输入下的方向接口一致性，`validation/scope/` 保留 SCOPE 对比入口。共用代码的方向一致性不等同于独立的温度/通量精度验证；迭代残差小也只证明所实现方程的数值闭合。当前模型是单点时间序列和统计几何近似，不是空间栅格、DEM 温度场或完整水文模型。研究比较应保存实际输入快照、材料光谱、气象、输出和中间诊断。
