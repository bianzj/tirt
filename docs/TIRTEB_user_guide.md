# TiRT 能量平衡操作说明

本文对应当前源码实现。已有安装版应核对其随附源码和版本；本次源码物理修正不会自动进入以前生成的安装包。

## 安装版

安装后启动 `TiRT`，点击 `TIRT-Dynamic` 进入能量平衡时间序列页面。带有内置运行时和 `data/` 的安装版无需另装 Python；批量运行或修改物理参数时使用同版本源码与数据。
## 源代码运行

~~~bat
start_tirt.bat
~~~

或：

~~~bash
python gui/server.py --port 8765
~~~

浏览器打开 [时间序列页面](http://127.0.0.1:8765/time.html)。命令行运行单个案例：

~~~bash
python tirteb_run.py --input input_tirteb.csv --output cases/my_case/output.csv
~~~

单场景（一个观测方向、一个波段，逐时间步输出组分温度、通量与亮温）：

~~~bash
python -m base.single_point --vza 0 --vaa 0 --wavelength 10.5
~~~

## GUI 操作顺序

1. 在 Site and period 设置经纬度、DOY 范围和气象文件。
2. 在 Observation 设置方向和观测时间筛选。
3. 在 Scene structure 选择平面、单坡、复合地形或城市，以及裸土、均质、垄行或树冠。
4. 在 Spectral information 设置观测波段和各组分 ρ、τ，`-` 表示取模型或实测值。表格显示的是观测波段发射率；全波段短波和宽带长波能量参数见下文。再选择叶片/土壤模型光谱或实测材料光谱。
5. 在 Plant physiology 和 Thermal/energy balance 中调整方法和参数。
6. 点击 Run time series，等待计算完成，再检查 `output.csv` 的 `closed`、能量残差和稳定度残差；任务执行完成不代表每个气象步都物理闭合。
7. 使用图表查看 Tv_C、LE、H 或 G，保存 `output.csv`、`observations.csv`、`intermediate.csv` 和实际输入快照。

## 输入文件

`input_tirteb.csv` 中与 TiRT `input.csv` 含义相同的参数（地表、坡面、地形、建筑、垄行、树冠 `crown_rad_a/crown_rad_b/crown_height/crown_density`、波段参数表）使用同名键；TiRT 的 `slope_angle`、`slope_aspect` 作为 `terrain_slope`、`terrain_aspect` 的同义键接受。两者并非完全同一套参数：TIRTEB 另有气象、时间、能量平衡、生理和热参数；观测方向用 `direction_source`、`vza`、`vaa`、`directions_file`，不支持 TiRT 的 `geometry_mode`、`geometry_file`、`geometry_pair_mode`、`principal_vza*`（出现时报错）。旧键 `crown_radius`、`crown_width` 报错，需改用 `crown_rad_a`（水平半径）、`crown_rad_b`（垂直半轴）、`crown_height`（冠底离地高度，与旧版含义不同）。不认识或不起作用的参数名（拼写错误，或已删除的 `G`、`lidfa`、`lidb`、`aero_zo`、`aero_d`、`aero_rbc`、`aero_cr`、`aero_rwc`、`urban_street_radiation_factor`、`doy`、`soil_temperature_c`、`urban_temperature_c`、`urban_temperature_method`、`soil_depth_step`、`urban_depth_step`、`width`、`stand_density` 等；树冠用 `crown_rad_a`、`crown_density`，垄行用 `row_*`）读取时给出带键名的警告并被忽略。参数说明见 [TIRTEB 理论与方法说明](TIRTEB_theory_and_methods.md#参数体系)。

`observation_time_mode=hour` 时，`observation_hours` 填当地钟点（0–24 小时，可带小数，如 `12;13.5`），每天相同钟点输出；非法值报错。

运行时在结果目录写入实际输入快照 `input.csv`（首行为 `# Effective TIRTEB input` 标记）。若该位置正是输入文件本身，或已有一个不是本程序生成的 `input.csv`，快照改写为 `input_used.csv`，不会修改输入文件。

方向文件（`direction_source=file` 时的 `directions_file`）每行至少包含：

~~~text
vza vaa
~~~

波段参数表文件（`spectral_source=file`）与 TiRT 格式相同，共 7 行；TIRTEB 中 `-` 表示该格取模型或实测值，例如：

~~~text
波段       0.55;0.86;8.6;10.5
叶片反射率 -;0.45;0.02;0.015
叶片透过率 -;0.40;0;0
土壤反射率 -;-;0.08;0.045
街道反射率 -;-;0.06;0.045
墙壁反射率 -;-;0.08;0.07
屋顶反射率 -;-;0.07;0.05
~~~

叶片可见光近红外文件（`vegetation_spectrum_source=measured`）可使用三列 wavelength、reflectance、transmittance；土壤文件（`soil_spectrum_source=measured`）使用两列 wavelength、reflectance；`urban_spectrum_file` 为七列（波长、屋顶/墙面/街道可见光近红外反射率、屋顶/墙面/街道发射率），其发射率在观测表对应格为 `-` 时用于观测，并在固定 10.5 µm 处提供城市宽带灰体默认参考；显式宽带覆盖另见下文。文件波长可以用 µm，也可以使用 nm；程序会自动识别小于等于 20 的波长为 µm。

## 气象与面积单位

气象文件第一行为记录数，随后每行至少 7 列数值；以下仅说明列顺序，文件中不写这一文字表头：

~~~text
time_DOY  air_temperature_C  vapor_pressure_hPa  pressure_hPa  wind_m_s  shortwave_W_m2  longwave_W_m2
~~~

气象时间按输入的当地时间和 `utc_offset_hours` 计算太阳位置，必须严格递增。所有气象值都按该时间戳时刻的瞬时值使用；原始数据是时段均值时也直接当作瞬时值，模型不平移时间。首步采用 `dtime` 秒（默认 1800）；后续步长取实际气象时间差，不要求等间隔。

自带的 `data/defined/meteo.txt` 是黑河流域张掖附近的观测，模板站点默认为大满站（`latitude=38.8555`、`longitude=100.3722`、`dem=1556`，北京时间 UTC+8）。经纬度必须与气象数据所在站点一致：经度写错时，太阳位置会与观测短波错开（例如误用北京经度，太阳正午会提前约 1 小时，傍晚晴空指数 kt 会大于 1）。`dem` 只作记录，空气密度使用气象文件中的气压计算。

气象短波按水平面总辐照度输入。单坡的方向短波反射率（`green`、`nir` 和 observations.csv 的短波列）也以水平面总入射为分母：坡面局地直射与天空反射率分别乘 `(1−f)·cos i/cos sza` 与 `f·(1+cos β)/2` 后相加，与 TiRT 相同。`shortwave_diffuse_fraction=0.2` 默认直射/漫射为 0.8/0.2，设为 −1 时采用 Erbs 估计，0–1 为显式漫射比例。单坡会自动执行水平直射到坡面的投影和天空视域修正，不要重复输入已经投影的坡面短波。

| 场景 | Rn、LE、H、G 的单位面积 | LAI 的基准 |
| --- | --- | --- |
| plane | 单位平面地面 | 同一地面 |
| slope / terrain_type=1 | 单位局地坡面 | 坡面面积 |
| terrain / terrain_type=2 | 单位水平地面 | 水平面积，内部换算到面群体 |
| urban | 单位水平地面，保留真实墙面积 | 单位街道冠层地面 |

城市叶片面积与完整街道面积同时参与能量预算；天底叶片可视比例不能作为街道能量面积的排除比例。复合地形使用 12 个圆锥方位面群体和 4 个共享组分温度，包含投影、统计遮阴和面间交换；它不是每个坡面具有独立温度的 DEM 模型。

## 宽带长波与观测波段

7 行波段参数表只控制对应观测波段。增加波段、删除波段或局部加密波段不会改变宽带能量平衡。默认宽带参数是固定 10.5 µm 的灰体参考，而非观测热红外波段平均。下列可选标量可以添加到 CSV，独立指定宽带参数：

~~~csv
thermal,broadband_leaf_emissivity,0.985,gray leaf emissivity
thermal,broadband_leaf_transmittance,0,gray leaf transmittance
thermal,broadband_soil_emissivity,0.955,gray soil emissivity
thermal,broadband_urban_emissivity_roof,0.95,gray roof emissivity
thermal,broadband_urban_emissivity_wall,0.93,gray wall emissivity
thermal,broadband_urban_emissivity_street,0.955,gray street emissivity
~~~

发射率和透过率必须为有限的 0–1 值，叶片二者之和不超过 1。上述示例会固定宽带值；若希望使用实测城市 10.5 µm 或土壤含水量发射率的默认参考，删除相应宽带覆盖行。逐波段的反射率、透过率仍由观测表保存和使用。

## 生理、含水量和导热

`stomatal_method=farquhar` 或 `ballberry`。两类叶片使用各自温度和 PAR：光照叶接受直接和对应散射 PAR，阴影叶仅接受散射 PAR。两种方案的气孔导度都不低于 `stomatal_g0`；Farquhar 在净光合 ≤0 时取该值，不再完全关闭气孔。`stressfactor` 是给定的植被胁迫参数。每段气动阻抗的上限由 `max_aerodynamic_resistance` 设置（默认 400 s m⁻¹，与 SCOPE 相同；0 表示不限制）。

`soil_moisture` 和 `satwater` 统一为体积含水量 m³ m⁻³（0.25 为 BSM 的容量参考值）；`soil_moisture` 读入后限制在 0.05–0.55，`satwater` 为正数，默认 satwater=0.45。为兼容旧文件，任一土壤含水量输入大于 1 时按百分数除以 100，模型内部只使用 m³ m⁻³，不再出现百分数形式。未显式给定 `soil_rss`/`rss` 时，默认阻抗为 `exp(8.206−4.255·soil_moisture/satwater)`。模板中已有 `soil_rss=200` 这种显式值时，默认保留常数含义。要明确改为含水量关系，可添加：

~~~csv
soil_thermal,soil_rss_method,sellers,moisture dependent surface resistance
~~~

若需要固定值，设 `soil_rss_method=constant` 并填写 `soil_rss`。这一关系不引入降水、渗流和根系吸水的水量演变；输入水分与植被胁迫仍为给定值。当前城市屋顶、墙面、街道的 LE 为零，城市潜热来自叶片，不能用于湿街道蒸发的直接模拟。

热参数单位必须区分：`soil_cs` 为质量比热 J kg⁻¹ K⁻¹，乘 `soil_rhos` 后才是体积热容；`urban_cs` 已是 J m⁻³ K⁻¹。城市街道的 G 和温度剖面推进均使用土壤热参数，屋顶/墙面使用城市热参数。导热支持 `direct`（净辐射比例）、`surface`（真实时间间隔的 7 小时热历史）和 `depth`（向后欧拉内部温度剖面）。depth 使用当前表面边界，不截断 Fourier 数；固定底边界或固定室内温度是模型假定。

## output.csv 列说明

每个气象时间步一行。数值为 `nan` 表示该量不存在或无定义（例如没有叶片时的叶温、夜间的短波反射率、非城市场景的屋顶温度）。每个地物都保留光照和阴影两个组分温度；某组分当前面积权重为 0（如夜间的光照组分）时温度照常输出，但不计入场景通量，可由对应 `frac_*_sunlit` 判断。

| 列名 | 单位 | 含义 |
| --- | --- | --- |
| `node` | – | 时间步序号（从 0 开始） |
| `time` | DOY | 气象时间（当地标准时间的年积日小数） |
| `veg` | – | 实际模型：`bare`、`homogeneous`、`row`、`crown` 或 `urban` |
| `surface` | – | `surface_model`：`plane`、`terrain` 或 `urban` |
| `Ta_C` | °C | 气温 |
| `Tss_C`、`Tsh_C` | °C | 光照、阴影土壤温度；城市中为光照、阴影街道温度（与 `Tstreet_*` 相同） |
| `Tls_C`、`Tlh_C` | °C | 光照、阴影叶片温度；LAI=0（裸土、城市 bare）时为 `nan` |
| `Rn`、`LE`、`H`、`G` | W m⁻² | 场景净辐射、潜热、感热、传导热（面积基准见上文表格）；G 为土壤（城市为屋顶、墙面、街道）传导热，正值表示进入材料 |
| `Troof_sunlit_C`、`Troof_shaded_C` | °C | 光照、阴影屋顶温度；仅城市 |
| `Twall_sunlit_C`、`Twall_shaded_C` | °C | 光照、阴影墙面温度；仅城市 |
| `Tstreet_sunlit_C`、`Tstreet_shaded_C` | °C | 光照、阴影街道温度；仅城市 |
| `Tleaf_sunlit_C`、`Tleaf_shaded_C` | °C | 城市街道冠层光照、阴影叶片温度；仅城市且 LAI>0 |
| `frac_soil_sunlit` | – | 土壤（城市为街道）能量面积权重中光照组分所占比例，阴影比例 = 1 − 该值 |
| `frac_leaf_sunlit` | – | 叶片光照比例 `LAIsun/LAI`（城市另计建筑遮阴）；LAI=0 时为 `nan` |
| `frac_roof_sunlit`、`frac_wall_sunlit`、`frac_street_sunlit` | – | 屋顶、墙面、街道面积中光照部分比例；非城市为 `nan` |
| `Tv_C` | °C | 天底方向 10.5 µm 亮温（TiRT 方向模型，使用该时刻组分温度） |
| `green` | – | 天底方向 0.55 µm（绿光）反射率；无太阳直射或短波为 0 时为 `nan` |
| `nir` | – | 天底方向 0.86 µm 近红外反射率；无日照时为 `nan` |
| `closed` | – | 组分能量与稳定度是否联合闭合 |
| `iterations` | – | 累计标量温度迭代数 |
| `energy_balance_residual_W_m2` | W m⁻² | 最大组分能量残差 |
| `L` | m | Obukhov 长度（中性可为 inf） |
| `ustar` | m s⁻¹ | 摩擦速度 |
| `stability_closed`、`stability_iterations`、`stability_residual_m_inv` | –、–、m⁻¹ | 稳定度外层是否闭合、外层迭代数、`1/L` 残差 |
| `stability_limited` | 0/1 | 1 表示稳定参数达到上限、使用了受限的有效 `1/L` |
| `stability_flux_area_factor` | – | 稳定度计算把场景通量换算到单位水平面积的系数（单坡为坡面/水平面积比，其余为 1） |
| `fallback_solver_used` | 0/1 | 1 表示该步启用了有界信赖域后备求解 |
| `temperature_solver` | – | 温度求解器名称；未启用后备求解时为 `relaxation`（固定松弛迭代） |
| `clearness_index` | – | 晴空指数 kt = 气象总短波 /（1361·E0·cos sza）；太阳在地平线下为 `nan` |
| `clearness_index_exceeds_one` | 0/1 | 1 表示 kt>1。气象值（含时段均值）一律按时间戳的瞬时值使用，日出、日落附近出现属预期；运行结束时另有汇总提示 |
| `solar_direct_normal_limit_W_m2` | W m⁻² | 直射法向辐照度上限（Meinel 晴空 DNI） |
| `shortwave_direct_normal_W_m2` | W m⁻² | 实际采用的直射法向辐照度 |
| `air_density_kg_m3` | kg m⁻³ | 空气密度 |
| `aerodynamic_roughness_length`、`aerodynamic_displacement` | m | 动量粗糙度、零平面位移；城市为 `nan` |
| `aerodynamic_canopy_weight` | – | 冠层与裸土空气动力学参数的混合权重（由 LAI 决定，0 为裸土、1 为完整冠层）；城市为 `nan` |

闭合相关字段的检查方法见下节。多方向、多波段观测写入 `observations.csv`（列 `node,time,veg,surface,vza,vaa,sza,saa` 加每个波段一列，列名为波长标签，如 `0p67`、`10p5`；短波为反射率，热红外默认亮温 K）。

## 联合闭合与输出检查

默认 `max_iter=50` 是每个外层稳定度试点的一批内层温度迭代数，`max_stability_iter=60` 为外层上限；累计 `iterations` 可以超过 50。温度使用固定松弛系数 0.08，单次变化限制为 ±1.5 K，不是自适应松弛或完整耦合雅可比求解。

外层求解 `1/L`，在当前温度能量预算充分闭合后，根据面积积分感热 H 和摩擦速度 ustar 更新 L；有符号区间时采用有保护的割线法，没有区间时采用当前值与目标值的半步更新。H>0 对应 L<0，H<0 对应 L>0，中性允许 L=inf。对流速度加入有效风速，基础风速下限由 `min_wind_speed` 设置（默认 0.05 m s⁻¹，作用于气象原始风速）。

| 输出字段 | 含义与检查 |
| --- | --- |
| `closed` | 组件能量和稳定度联合闭合；应逐步检查 |
| `energy_balance_residual_W_m2` | 最大组件能量残差；最终公共阈值为 5 W m⁻²，稳定度采样使用更严格的内层 1e−4 阈值 |
| `L` | 实际返回的 Obukhov 长度 m，中性可为 inf |
| `ustar` | 摩擦速度 m s⁻¹ |
| `stability_closed` | `abs(1/L−target) ≤ 1e−7+1e−6·abs(target)+|c|·1e−4/ustar³`（最后一项为 1e−4 W m⁻² 的 H 容差传递到 1/L 的量，c 为浮力系数） |
| `stability_residual_m_inv` | 上述倒数长度残差，单位 m⁻¹ |
| `stability_iterations` | 外层稳定度试点数 |
| `iterations` | 累计标量温度迭代及有界信赖域后备求值数 |

启用 `output_intermediate=yes` 后，`intermediate.csv` 还包含组分 `EB_residual_*`、阻抗、能量面积权重、`stability_zeta`、`stability_limited`、`L_target` 和城市/地形诊断。稳定函数先按参考高度统一限制整个廓线使用的有效 `1/L`，使参考高度稳定参数不超过 1；各高度共用这一有效倒数长度，不分别裁剪两端订正。实际 L 仍由 H 求出且不裁剪。`stability_limited=True` 表示进入限制范围，不能据此声称强稳定夜间交换已独立验证。城市 `urban_SW_absorbed + urban_SW_escaped` 应等于其外部短波入射。

`closed=False` 的时间步仍可输出并进入后续时间步；应保留记录并检查具体残差、边界条件和参数。数值闭合、辐射守恒、与 TiRT 的方向接口一致性和独立实测精度属于不同层次，不能互相替代。详细方程和近似见 [理论与方法说明](TIRTEB_theory_and_methods.md)。

## 绘图

安装 matplotlib 后：

~~~bash
python analysis/plot_observations.py --result-dir cases/my_case --time 214.5 --bands 0p67 0p86 10p5
~~~

图件写入案例目录的 `plots/`。0p67、0p86 和 10p5 是结果列名中的波长标签。

## 常见问题

- 没有观测记录：检查 observation_time_mode 和气象文件中的实际时间是否匹配。
- 叶片或土壤光谱文件报错：至少提供两行有效数值，并保证波长列和光谱列顺序正确。
- 复合地形报错 “Total coverage terrain_density * pi * terrain_radius^2 must be less than 1”：`terrain_density` 是每平方米的圆锥数，需满足 `terrain_density·π·terrain_radius² < 1`。
- 计算完成但 `closed=False`：先区分最大组件能量残差与 `stability_residual_m_inv`，再检查实际步长、热参数单位、强稳定限制、辐射入射和温度保护；不要仅因任务状态 Complete 就认定闭合。
- 波段修改导致亮温变化：逐波段光学性质会改变观测，宽带能量参数独立；对比温度时应同时保存宽带覆盖行和土壤含水量设置。
- 单坡/复合地形比较：先统一局地坡面或水平地面面积基准，再比较通量；复合地形仅有四个共享温度。
