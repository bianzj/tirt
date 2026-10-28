# TiRT / TIRTEB 理论与方法说明

本说明给出当前标准版（`canopy_closure=standard`、TIRTEB `stomatal_method=farquhar`）默认路径上实际执行的公式，包含 TiRT 方向辐射模型与 TIRTEB 能量平衡模型两部分。每个公式给出符号定义、实现位置和默认参数；实验性选项（`legacy`、`upgrade` 及旧几何/光学分支）不在本说明范围内。

## 1 总体结构与调用路径

TiRT 计算给定组分温度下的方向辐亮度/亮温与方向反射率；TIRTEB 在 TiRT 之上求解能量平衡、组分温度、阻抗与通量的时间序列：

- TIRTEB 主链：气象 → 短波/长波吸收 → 能量平衡 → 组分温度 → 通量；
- TIRTEB 方向观测：组分温度与观测方向交给 TiRT，得到方向辐亮度。

$$
L_{obs} = \text{TiRT}\left(T_{ss}, T_{sh}, T_{ls}, T_{lh};\ \lambda,\ \theta_v,\ \phi_v\right)
$$

TiRT 的调用路径：

```text
run.run → _scenario_configs → _simulate_rows
  λ < 2.5 µm → _optical_band
       slope        → optical.reflectance：局地平面递归（直射项 + 天空项）
       urban        → brightness.urban_directional_reflectance
       plane  hom/row → canopy_closure.reflectance
       plane  crown   → optical_crown.crown_reflectance → canopy_closure.reflectance
       terrain 植被   → optical_terrain_standard.reflectance
       terrain 裸地   → vnir.canopy_reflectance（LAI=0，见 §13）
  λ ≥ 2.5 µm → _surface_result
       plane  → canopy_closure.thermal_weights
       slope  → 局地几何下的 plane
       terrain → TerrainForest（裸地为 Terrain）+ 节点散射
       urban  → Urban / Urban_Veg → urban_radiosity
```

TIRTEB 的调用路径：

```text
tirteb_run.run → 逐步 meteo → pipeline.run_node
  太阳几何 solar.solar_angles
  短波分配 tirteb_radiation.prepare_shortwave（Erbs、Meinel 晴空上限、晴空指数）
  稳定度外层迭代（x = 1/L）
     → 气动阻抗 aerodynamic.aerodynamic_resistance
     → 组分吸收 tirteb_radiation（四组分 × 节点吸收）
     → 生理 biochemical.suresist（farquhar）
     → 能量平衡 eb.rebalance（各组分温度）
  → 场景通量 scene_fluxes → 热历史 evapo.update_thermal_history → 输出
```

## 2 通用定义

**方向与角度**。太阳天顶角 θs、方位 φs；观测天顶角 θv、方位 φv；μ = cos θ。相对方位折叠到 [0°, 180°]：

$$
\phi = \left| (( \varphi_v - \varphi_s + 180 ) \operatorname{mod} 360) - 180 \right|
$$

**叶投影系数**。球形叶倾角分布取 G = 0.5，方向消光系数：

$$
k(\theta) = \frac{G}{\mu} = \frac{0.5}{\cos\theta}
$$

**方向间隙（Beer–Lambert）**：

$$
P(\theta) = \exp\left(-k(\theta)\,L\right) = \exp\left(-\frac{0.5\,L}{\cos\theta}\right)
$$

**太阳—观测角距离**（热点相关）：

$$
\Delta = \sqrt{(\tan\theta_v - \tan\theta_s)^2 + 4\tan\theta_v\tan\theta_s\sin^2(\phi/2)}
$$

**叶尺度热点相关**。`hspot` 为叶片尺度与冠层厚度之比：

$$
a = \frac{2\Delta}{h_{spot}\,(k_s + k_v)},\ w = \frac{1 - e^{-a}}{a}
$$

Δ = 0 时 w = 1；`hspot` ≤ 0 时关闭相关（w = 0）。w 用于把两条射线的共同覆盖深度从独立值收缩：

$$
t_{sh} = \min\left(t_v,\ t_s,\ \sqrt{t_v t_s}\,w\right)
$$

## 3 太阳与观测几何

**太阳位置**采用 NOAA 太阳位置算法（含大气折射），输入为当地标准时、`utc_offset_hours`、经纬度和日序。天顶角用于辐射计算的截断：

$$
\theta_{v,model} = \min(|\theta_v|,\ 89.999^\circ),\ \theta_{s,model} = \min(|\theta_s|,\ 89.9^\circ)
$$

**坡面局地变换**（坡度 β、坡向 A，Δ = φ − A）：

$$
x,\ = \cos\Delta\sin\theta\cos\beta - \cos\theta\sin\beta
y,\ = \sin\theta\sin\Delta
z,\ = \sin\theta\cos\Delta\sin\beta + \cos\beta\cos\theta
\theta',\ = \arccos z,\ (\theta' > 90^\circ  180^\circ - \theta'),\ \phi' = \operatorname{atan2}(y, x)
$$

坡面入射余弦与可见判据：

$$
\cos i = \cos\beta\cos\theta + \sin\beta\sin\theta\cos(\varphi - A)
$$

**观测方向组合**（geometry_mode=0）：观测天顶角 × 观测方位 × 太阳天顶角 × 太阳方位 的交叉积；默认 7×4×1×1 = 28 个方向（θv 0–60°、步长 10°，φv 0/90/180/270°，θs 30°、φs 0°）。

## 4 四组分可见比例

四组分为日照土壤 P_ss、阴影土壤 P_sh、日照叶片 P_ls、阴影叶片 P_lh，满足

$$
P_{ss} + P_{sh} + P_{ls} + P_{lh} = 1
$$

三者共用同一构造：先求方向间隙 P_v、P_s 与联合间隙 P_ss，再由日照叶近似给出 P_ls，其余为阴影。

### 4.1 裸土

$$
\mathbf{P} = (1,\ 0,\ 0,\ 0),\ L = \varepsilon_s(\theta_v)\,B(\lambda, T_{soil,sun})
$$

### 4.2 均质冠层

均质冠层等价于行间空隙为零的垄行（w=0.5、b=0、h=1），此时两相替代退化为 c = 0、t = kL：

$$
P_v = e^{-k_v L},\ P_s = e^{-k_s L}
$$

$$
P_{ss} = \exp\left[-\left(t_s + t_v - \min(t_s,\ t_v,\ \sqrt{t_s t_v}\,w)\right)\right],\ P_{sh} = P_v - P_{ss}
$$

### 4.3 垄行冠层

垄行由行宽 w、行间空地 b、行高 h、行向 φ_r 描述，行内叶密度

$$
\rho_L = \frac{L\,(w+b)}{w\,h}
$$

行向间隙（方位与行向平行或天底观测）：

$$
P = \frac{w\,e^{-G\rho_L h/\mu} + b}{w + b}
$$

一般方向把行结构投影为一维周期条带：w' = w/|sin φ_r|、b' = b/|sin φ_r|、h' = h tan θ、n = ⌊h'/(w'+b')⌋，按条带交叠情形逐段积分，单段核为

$$
S(c_1,x_1,c_2,x_2,k) = c_1 e^{-x_1} + c_2 e^{-x_2} - k\,e^{-x_2}\left(e^{x_2-x_1} - 1\right)
$$

求得的 P 再作两相替代，把冠层拆成"覆盖相"与"空相"：

$$
(1-c)\,t = kL,\ c + (1-c)\,e^{-t} = P
$$

两条射线在行间的错位与覆盖条带重叠：

$$
s = h\tan\theta\sin\varphi_r
O = \frac{\operatorname{overlap}\left((1-c_s)P_r,\ (1-c_v)P_r,\ (s_s-s_v)/2,\ P_r\right)}{P_r}
$$

$$
P_{ss} = (1-a-b+O) + (a-O)e^{-t_s} + (b-O)e^{-t_v} + O\exp\left[-(t_s+t_v - t_{sh})\right]
$$

其中 a = 1 − c_s、b = 1 − c_v，P_r = w + b，最后取 P_ss = min(P_ss, P_s, P_v)。该闭式与截面数值积分一致（差 < 1e-5）。

### 4.4 离散树冠

树冠为泊松分布的浑浊旋转椭球：株密度 n、水平半径 R、垂直半轴 H、冠底高 h_b，冠内消光系数

$$
k = \frac{G\,L}{n\,V},\ V = \frac{4}{3}\pi R^2 H
$$

单冠投影面积与平均光学厚度：

$$
A(\theta) = \pi R\sqrt{R^2 + (H\tan\theta)^2},\ \tau(\theta) = \frac{k\,V}{A\cos\theta}
$$

冠体平均透过率（x = 1.5τ）：

$$
\bar t = \frac{2\left[1-(1+x)e^{-x}\right]}{x^2}
$$

单冠遮挡面积 J = A(1 − t̄)，土壤可见比例

$$
P_v = e^{-n J_v}
$$

**冠尺度热点**由两投影透镜重叠给出。冠心高度 H_c = max(h_b,0) + H，两投影中心距

$$
d = H_c\left|\tan\theta_v\,\hat e_v - \tan\theta_s\,\hat e_s\right|
$$

冠尺度相关权与联合透过率：

$$
a_c = \frac{2\Delta}{(k_s+k_v)}\cdot\frac{4H/3}{2H\,h_{spot}},\ w_c = \frac{1-e^{-a_c}}{a_c}
$$

$$
t_{vs} = \min\left\{\exp\left[\min\left(\ln t_v + \ln t_s + \min(\sqrt{\tau_v\tau_s}\,w_c,\ \tau_v,\ \tau_s),\ 0\right)\right],\ t_v,\ t_s\right\}
$$

共同遮挡（两圆透镜重叠面积 × 一冠同时挡住两射线的概率）：

$$
K = \operatorname{Lens}\left(\sqrt{A_v/\pi},\ \sqrt{A_s/\pi},\ d\right)\cdot\max\left(1 - t_v - t_s + t_{vs},\ 0\right),\ K \leq \min(J_v, J_s)
$$

$$
P_{ss} = P_v\exp\left[-n\,(J_s - K)\right],\ P_s = \left. P_{ss} \right|_{\theta_v = \theta_s}
$$

注意冠尺度热点的深度尺度为 (4/3)H，而 §4.5 的叶尺度相关用 2H，两者在同一 `hspot` 下相差 1.5 倍（见 §13）。

### 4.5 覆盖相双层日照叶近似（均质/垄行/树冠共用）

日照叶比例由"覆盖相内的双层近似"给出。记 P_v 的两相参数为 (c_v, t_v)、P_s 为 (c_s, t_s)，叶尺度共享深度 t_sh 见 §2：

$$
L_c = \frac{L}{1-c_v},\ p_v = e^{-t_v},\ p_s = e^{-t_s},\ p_{ss} = e^{-(t_v+t_s-t_{sh})}
$$

在覆盖相内（叶面积 L_c）作双层近似：上层厚度比例 u，下层 1−u，

$$
k_v = -\frac{\ln p_v}{L_c},\ k_s = -\frac{\ln p_s}{L_c},\ o = \operatorname{clip}\left(\frac{\ln\left(p_{ss}/(p_v p_s)\right)}{L_c},\ 0,\ \min(k_v,k_s)\right)
$$

$$
\iota = \max(0.82 - 0.08\,L_c,\ 0.02)
u = \sqrt{\left[\frac{-\ln\left(1-(1-p_v)\iota\right)}{k_v L_c}\right]\left[\frac{-\ln\left(1-(1-p_s)\iota\right)}{k_s L_c}\right]}
$$

$$
P_{ls}' = 1 - e^{-k_v L_c u} + \left(1 - e^{-o L_c(1-u)}\right)e^{-(k_s+k_v-o)L_c u}
$$

$$
P_{ls} = \operatorname{clip}\left((1-c_v)\,P_{ls}',\ 0,\ 1-P_v\right),\ P_{lh} = 1 - P_v - P_{ls}
$$

其中 ι 为经验截获比（"上层可见叶全部受照"的双层近似），是当前模型中唯一按经验系数构造的叶尺度几何；该式对太阳与观测不对称（§13）。

## 5 热红外辐射与发射率

**基尔霍夫发射率**由叶片与土壤的反射率、透过率给出：

$$
\varepsilon_l = \max(1-\rho_l-\tau_l,\ 0)
\varepsilon_s = 1-\rho_s
\varepsilon_{roof} = 1-\rho_{roof},\ \varepsilon_{wall} = 1-\rho_{wall},\ \varepsilon_{street} = 1-\rho_{street}
$$

**普朗克辐亮度**（W m⁻² sr⁻¹ µm⁻¹，λ 单位 µm）：

$$
B(\lambda,T) = \frac{c_1\times10^{4}}{\lambda^5\left[\exp\left(\dfrac{c_2}{\lambda T}\right)-1\right]},\ c_1 = 11910.439,\ c_2 = 14388.291\ \text{µm K}
$$

**反普朗克（亮温）**：

$$
T_b = \frac{c_2}{\lambda\ln\left[\dfrac{c_1\times10^{4}}{L\,\lambda^{5}} + 1\right]}
$$

**方向辐亮度**为各组分有效发射权重乘以组分温度的普朗克辐亮度：

$$
L(\lambda,\theta_v) = \sum_{i\in\{ss,sh,ls,lh\}} w_i(\theta_v)\,B(\lambda,T_i)
$$

平地冠层的有效权重为直接可见部分加节点散射贡献：

$$
w_i = P_i\,\varepsilon_i + S_i
$$

其中 S_i 为 §6 的节点多次散射权重（只计二阶及以上）。**组分温度**由输入给出（`*_c` 摄氏度或 `*_k` 开尔文）。

**天空反射项**在给定 `sky_temperature_k` 时叠加（默认黑天空，不执行）：

$$
L' = L + \max(0,\ 1-\varepsilon_{eff})\,B(\lambda,T_{sky})
\varepsilon_{eff} = \frac{L_{300}}{B(\lambda, 300\,\mathrm{K})}
$$

## 6 节点辐射度多次散射

### 6.1 共同框架

把冠层与土壤离散为若干节点，节点容量 A_i 为土壤面积比或叶面积 4G L_LAI。节点间交换因子满足互易

$$
A_i F_{ij} = A_j F_{ji},\ e_i = \operatorname{clip}\left(1 - \sum_j F_{ij},\ 0,\ 1\right)
$$

扩散截获用指数积分 E₃：

$$
X(\tau) = 1 - 2E_3(\tau)
$$

每个节点有一个上行态；冠层节点另有一个下行态（发射各半）。从入射到出射的散射矩阵 S，叶片按球形双朗伯：

$$
b = \frac{2\rho+\tau}{3}\ (),\ f = \frac{\rho+2\tau}{3}\ ()
$$

全阶解由线性方程给出（J₁ 为一次散射后的出射功率）：

$$
\mathbf{J} = \left( \mathbf{I} - \mathbf{S} \mathbf{A}^{\mathrm{T}} \right)^{-1} \mathbf{J}_1
$$

传感器权重为节点的可见比例乘其上行辐亮度，半球积分用 32 点 Gauss（μ）× 24 方位（一个象限，行结构对称）。

**热红外散射**：节点日照比为土壤 sun_i/A_i、叶片 sun_i/(k_s L_i)，源为 ε·σ 与 ε(1−σ)；散射权重取

$$
S_i = W_i\left[\mathbf{J}(\rho_s = 1-\varepsilon_s,\ \rho_l = 1-\varepsilon_l,\ \tau_l = 0) - \mathbf{J}_0\right]
$$

### 6.2 均质冠层 hom3（土壤 + 2 个等厚冠层）

层厚 τ = kL/2，各阶交换因子：

$$
A F_{ii} = 4\tau - 2X(\tau)
A F_{i,i\pm1} = 2X(\tau) - X(2\tau)
$$

$$
A F_{ij} = \left[X((g+1)\tau) - X(g\tau)\right] - \left[X((g+2)\tau) - X((g+1)\tau)\right],\ g = |i-j|-1
$$

$$
A F_{0i} = X(i\tau) - X((i-1)\tau)
$$

方向分配为土壤 e^{−kL}，第 j 层 exp[−kL·(上方层数)/2]·(1 − e^{−kL/2})。

### 6.3 垄行冠层 row3（覆盖土壤 + 开放土壤 + 1 个冠层）

开放土壤比例 f_o = b/(w+b)。对每个方向用两相参数求分土壤半球间隙（share = min(1, c/f_o)）：

$$
M_{cov} = 2\!\int_0^1 \mu\,e^{-t}\,d\mu
M_{open} = 2\!\int_0^1 \mu\left[\text{share} + (1-\text{share})e^{-t}\right]d\mu
$$

$$
F_{soil_k\to canopy} = 1-M_k
F_{canopy\to soil_k} = \frac{A_k(1-M_k)}{4GL}
F_{cc} = 1 - \sum_k \frac{A_k(1-M_k)}{2GL}
$$

### 6.4 树冠 crown2（单一土壤 + 1 个冠层）

方向间隙用 §4.4 的冠体间隙（`hspot` 取 0.1，与方位无关），半球间隙 M = 2∫μ P(θ)dμ：

$$
F_{soil\to canopy} = 1-M,\ F_{canopy\to soil} = \frac{1-M}{4GL},\ F_{cc} = 1-\frac{1-M}{2GL}
$$

crown2 不区分覆盖与开放土壤，且只有一层冠层。

### 6.5 短波多次散射

**直射**（二阶及以上，单位水平直射按 split(θs) 进入各节点）：

$$
R_{beam} = \pi \sum_i W_i\left(J_i - J_{1,i}\right)
$$

**各向同性天空**（全阶，进入节点按逃逸权 e_i A_i 分配）：

$$
R_{sky} = \pi \sum_i W_i J_i
$$

合成（f_d 为散射比例）：

$$
R_{multi} = f_d R_{sky} + (1-f_d)\left[\theta_s<90^\circ\right]R_{beam}
$$

## 7 短波一次散射与方向反射率

**散射比例** f_d = `shortwave_diffuse_fraction`（默认 0.2），直射权重 D = (1−f_d)·[受照]。

**球形叶双向散射系数**（SAIL/Verhoef，64 点 Gauss 对叶倾角积分）：

$$
sob = \sum_j w_j\frac{\pi f_\rho(\theta_s,\theta_v,\phi,\theta_{l,j})}{\mu_s\mu_v}
sof = \sum_j w_j\frac{\pi f_\tau(\theta_s,\theta_v,\phi,\theta_{l,j})}{\mu_s\mu_v}
$$

其中 f_ρ、f_τ 为单叶的反射/透过散射相函数，θ_l 为叶倾角（μ 在 [0,1] 均匀 = 球形叶倾角分布）。

**方向反射率（BRF/HDRF）**，一次散射用 §4 的四组分，多次散射用 §6.5：

$$
R_1 = D\left[\frac{sob\,\rho_l + sof\,\tau_l}{k_v}P_{ls} + r_s P_{ss}\right]
$$

$$
R = R_1 + R_{multi}
$$

裸地或 L = 0 时退化为 R = D r_s + f_d r_s；坡面背光时 P → (0, P_ss+P_sh, 0, P_ls+P_lh) 且 D = 0。

节点方案：均质 hom3、垄行 row3、树冠 crown2（与热红外共用同一套节点）。

## 8 单坡

坡面先把全局方向变换到局地坐标（§3），再在局地几何下用 §4–§7 的平地公式；LAI 与株密度按坡面面积计（树木垂直于坡面生长）。

**热红外**：坡面背光（cos i < 0）时日照组分改用阴影温度；相对坡面法线天顶角大于 89° 的方向省略输出。

**短波**：以水平面总入射归一（与坡面能量平衡一致），只计各向同性天空视因子、不计周围地面反射：

$$
R = W_{dir}\,R_{local}(f_d=0) + W_{sky}\,R_{local}(f_d=1)
$$

$$
W_{dir} = (1-f_d)\frac{\max(\cos i,\ 0)}{\cos\theta_s},\ W_{sky} = f_d\frac{1+\cos\beta}{2}
$$

局地 R 用局地观测/太阳角与局地行向计算，局地太阳天顶角截断到 89.9°。

## 9 复合地形

地形为高度 h、半径 r、密度 n 的圆锥群体。锥底占地与锥影超出底面的部分（Li–Strahler）：

$$
p_{occ} = \sum n\pi r^2,\ L = h\tan\theta\ ( r),\ \gamma = \arcsin(r/L),\ \Lambda = n\left(\cot\gamma + \gamma - \frac{\pi}{2}\right)r^2
$$

地面可见与受照比例（Ω 为两影重叠比）：

$$
P_{g,v} = e^{-\Lambda_v}(1-p_{occ})
P_{g,vs} = e^{-(\Lambda_v + \Lambda_s\Omega)}(1-p_{occ})
\Omega = \frac{\left|\tan\theta_v\hat e_v - \tan\theta_s\hat e_s\right|}{\tan\theta_v + \tan\theta_s}
$$

锥面可见投影与高度带划分：

$$
A_{mount}(\gamma) = \pi r^2,\ L \leq r
\left(\cot\gamma + \gamma + \frac{\pi}{2}\right)r^2,\ L > r
$$

锥面分 3 个高度带，带权 r₁πΔh + 2r₁ + Δh，r₁ = (1−z/h)r；他山遮挡用 exp(−Λ_s Ω)。面元自照比例由 10 个方位面元按 max(n·v,0) 加权、以 n·s > 0 判定受照。

面元局地角与冠层：面元法线天顶 β₁ = π/2 − atan(r/h)，局地观测/太阳天顶角为 arccos(n·v)、arccos(n·s)，面元冠层用 LAI·cosβ₁ 与株密度·cosβ₁ 计算 §4–§6；背光面元全部按阴影处理。

**地形单次交换**（默认 `forest_transport=legacy`）：观测射线落在锥上的概率 i0v = 1 − exp(−Σ n A_mount(θ) r²/...)，面元上/下半球视因子 (1±sin α)/2（α = atan(r/h)），三条路径

$$
x_{ff} = p\,i0v\,\rho_f,\ x_{fg} = w_s\,\rho_g,\ x_{gf} = e_d\,i0v\,\rho_f
$$

其中 ρ_f、ρ_g 为面元与地面的冠层漫反照率。反射分量按 x 与日照份额分配到土壤与叶片，并在 Σscatter ≤ 1 − Σdirect 处按比例截断。

**复合地形短波**：一次散射部分直接使用面元的四组分（平地部分加面元部分，面元部分乘入射 max(cos i_f,0)/μs、他山遮挡与自照比例，并按可见锥面归一），多次散射用面元的局地节点算子与天空视因子；地形单次交换项把反照率换成冠层漫反照率。

## 10 城市

城市由等高屋顶、四面墙与街道组成：p_roof = ΣnLW，墙投影遮挡深度

$$
\Lambda_v = \sum n H\tan\theta_v\left(L\left|\cos(\varphi_v-\varphi_L)\right| + W\left|\cos(\varphi_v-\varphi_W)\right|\right)
$$

联合遮挡 Λ_vs = max(Λ_v,Λ_s) + min(Λ_v,Λ_s)·Ω。街道可见与受照：

$$
P_{st,v} = e^{-\Lambda_v}(1-p_{roof}),\ P_{st,sun} = e^{-\Lambda_{vs}}(1-p_{roof})
$$

屋顶（等高）全部可见且受照 P_roof,sun = p_roof；墙可见 P_wall = 1 − P_st,v − p_roof，墙受照按朝向与太阳、观测夹角判定并乘高度带平均的地形间隙。

**街道—四墙辐射度**（`urban_transport` 默认 radiosity）：节点为街道与四面墙，视因子由 20 点 μ × 72 方位的天空视因子构造，满足面积互易 A_iF_ij = A_jF_ji，天空份额为 1 − ΣF；求解

$$
\mathbf{J} = \left(\mathbf{I} - \operatorname{diag}(\rho)\mathbf{F}\right)^{-1}\mathbf{S}
$$

组分权重为直接可见部分加各节点的反射项；街道冠层按半球叶覆盖与叶日照比混合发射率与源项，冠层内部再用 §6 的节点散射。

**城市短波**：屋顶按 ρ_roof、街道按 ρ_street、墙按 ρ_wall 反射，墙间互遮用启发式因子，辐射度形式同热红外；有植被时街道边界的方向响应换成 §7 的冠层算子，节点反照率取冠层的漫反照率。

## 11 TIRTEB 能量平衡

### 11.1 输入与时间

气象文件每行一个时间戳（当地标准时，日序的小数形式）加 6 个气象量：气温、水汽压、气压、风速、入射短波、入射长波。所有气象值按该时间戳的瞬时值使用；原始数据为时段均值时也直接当作瞬时值。首步步长取 `dtime`（默认 1800 s），后续步长取实际气象时间差。

### 11.2 太阳几何与短波分配

太阳位置同 §3。晴空指数与散射分解：

$$
k_t = \frac{R_{sw}}{\max(S_0\cos\theta_s,\ \epsilon)},\ S_0 = 1367\,E_0
$$

其中 E_0 为日地距离修正。散射比例取请求值与晴空约束下的较大者：

$$
f_d = \max\left(f_{d,req},\ 1 - \min\left(\frac{DNI_{limit}\,\cos\theta_s}{R_{sw}},\ 1\right)\right)
$$

晴空法向直射上限（Meinel）与大气质量（Kasten–Young）：

$$
DNI_{limit} = 1367\,E_0\,0.7^{AM^{0.678}}
AM = \frac{1}{\cos\theta_s + 0.50572\,(96.07995^\circ - \theta_s)^{-1.6364}}
$$

直射超过上限的部分转为散射，总短波保持不变。

### 11.3 组分吸收

入射短波按散射比例分成直射与散射两支，再按太阳光谱权重 f_sun、天空光谱权重 f_sky 展布到 0.4–2.4 µm 的 2001 个 1 nm 波段：

$$
E^{sun}_{\lambda} = R_{sw}(1-f_d)f_{sun,\lambda}\cdot10^{-3},\ E^{sky}_{\lambda} = R_{sw} f_d f_{sky,\lambda}\cdot10^{-3}
$$

冠层被离散为与 TiRT 相同的节点（均质 hom3、垄行 row3、树冠 crown2）。设节点 i 的一次吸收为 J₁,i、多次吸收为 J_i^scat，则每个波段的节点吸收为

$$
A_{i,\lambda} = E^{sun}_{\lambda} J_{1,i} + E^{sun}_{\lambda} J^{scat}_{i} + E^{sky}_{\lambda}\left(J_{1,i}+J^{scat}_{i}\right)
$$

直射只落在受照组分上：按受照土壤面积、受照叶面积归一后得到单位受照面积的直射吸收；散射部分按土壤与叶片各自的总面积归一，再按受照/阴影面积比例分配，得到四组分（光照叶、阴影叶、光照土壤、阴影土壤）的吸收。四组分的净辐射为

$$
R_{n,i} = A^{direct}_i + A^{diffuse}_i + L_{n,i}
$$

其中 L_n 为长波项：长波在同一个互易节点几何上求全阶，节点发射按光照/阴影温度的四次方以面积混合；等温节点处于同温黑体辐射环境时净长波为零。

### 11.4 气动阻抗

风速先取下限再与对流速度合成：

$$
u_* = \frac{\kappa\,U}{\ln\frac{z-d}{z_0} - \psi_m(\zeta)}
U = \sqrt{\max(u,\ u_{min})^2 + w_*^2}
$$

其中 u 为气象文件原始风速（读入时不截断），u_min = `min_wind_speed`（默认 0.05 m/s）；自由对流速度取固定混合层高度 z_mix = 1000 m：

$$
w_* = \left(\frac{g\,z_{mix}\,H_{scene}}{\rho\,c_p\,T_a}\right)^{1/3},\ z_{mix} = 1000\ \text{m}
$$

只有向上的感热（H > 0）产生阵风，夜间 H ≤ 0 时 w* = 0。

稳定度函数（ζ = (z−d)/L）：

$$
\psi_m =
2\ln\frac{1+x}{2} + \ln\frac{1+x^2}{2} - 2\arctan x + \frac{\pi}{2},\ x=(1-16\zeta)^{1/4},\ \zeta<0
-5\zeta,\ \zeta>0
$$

冠层几何用 Raupach 关系，d 与 z0 由 LAI 与拖曳参数给出（与 SCOPE 相同）。各段的体阻抗为冠层上方段、叶边界层段与冠层内到土壤段之和，每段受上限约束：

$$
r_{aa} = \sum_{segments} \min\left(r_{segment},\ r_{max}\right),\ r_{max} = 400\ \mathrm{s\,m^{-1}}
$$

LAI < 0.5 时冠层与裸土的 u*、土壤阻抗按权重 w = s²(3−2s)（s = LAI/0.5）混合，使 LAI → 0 连续过渡到裸土。

### 11.5 稳定度外层迭代

以 x = 1/L 为未知量，用割线法配合包围区间求解

$$
x = \frac{\kappa g\,H_{scene}}{\rho c_p T_a\,u_*^3}
$$

每次给定 x 后重新求温度、气孔阻抗、气动阻抗、阵风与组分通量。收敛判据

$$
\left|x - x_{target}\right| \leq 10^{-7} + 10^{-6}\left|x_{target}\right| + \frac{|c| \cdot 10^{-4}}{u_*^{3}}
$$

最后一项把 H 的容差（1e−4 W m⁻²）传递到 1/L；外层内部用该容差的一半，使最终复核留有余量。内层温度不收敛时调用信赖域后备求解器，仍不满足则联合求解温度与稳定度。

### 11.6 能量平衡

每个组分（光照叶、阴影叶、光照土壤、阴影土壤）分别满足

$$
R_n - \lambda E - H - G = 0
$$

感热与潜热用体阻抗形式：

$$
H = \rho c_p \frac{T_c - T_a}{r_a}
\lambda E = \frac{\rho c_p}{\gamma}\cdot\frac{e_s(T_c) - e_a}{r_a + r_s}
$$

其中 ρ 为湿空气密度（Rd 与虚温 Tv = T/(1 − 0.378 e/p)），γ 为干湿表常数，e_s 用 Tetens 公式。组分温度由松弛迭代更新，

$$
T_c^{(k+1)} = T_c^{(k)} + w\frac{R_n - \lambda E - H - G}{\rho c_p/r_a + \rho c_p\,\Delta/(r_a+r_s) + 4\varepsilon\sigma T_c^{3}}
$$

松弛因子 w 随迭代次数递减；能量平衡残差以 5 W m⁻² 为判据（`RAD_THRESHOLD`）。

### 11.7 气孔阻抗

Farquhar 方案按光合速率求导度：Vcmax 随温度用 Arrhenius 型函数标定，电子传输用 Collatz 温度响应，C3/C4 分别取 q 与 kpep。净光合 A ≤ 0 时气孔导度取最小导度

$$
g_s = \max\left(\frac{1.6\,A\,p}{c_s - c_i},\ g_0\right)
r_{sw} = \frac{\rho}{M_{air}}\cdot\frac{1}{g_s}
$$

g₀ = `stomatal_g0`（默认 0.01 mol m⁻² s⁻¹），即不设完全关闭的气孔阻抗。土壤表面阻抗用 Sellers 型（默认）或常数：

$$
r_{ss} = \exp\left(8.206 - 4.255\,\frac{\theta_{soil}}{\theta_{sat}}\right)
$$

### 11.8 土壤与城市导热

模板默认 `soil_heat_method=surface`：把土壤视为半无限均匀介质，用地表温度历史按热惯性解析积分（Carslaw 解，回溯 48 h）：

$$
G(t) = \frac{2I}{\sqrt{\pi}}\int_{t-48\,\mathrm{h}}^{t}\frac{\mathrm{d}T_s/\mathrm{d}\tau}{\sqrt{t-\tau}}\,\mathrm{d}\tau,\ I = \sqrt{\rho_s c_s \lambda_s}
$$

历史按分段线性温度变化逐段积分，每段贡献 (ΔT/Δt)(√t₂ − √t₁)。

可选 `depth` 方法按一维热传导有限体积求解，深度 2 m（顶层层厚 0.01 m，按 1.3 倍逐层加厚），表面为 Dirichlet 温度、底部零通量：

$$
\rho_s c_s\frac{\partial T_s}{\partial t} = \frac{\partial}{\partial z}\left(\lambda_s\frac{\partial T_s}{\partial z}\right)
\left.\frac{\partial T_s}{\partial z}\right|_{z=2\,\mathrm{m}} = 0
$$

表面通量取顶层导热加半层蓄热。城市屋顶与墙面默认 `depth` 方法（`urban_heat_method=depth`），按材料厚度与室内固定温度 24 °C 求解围护结构导热；街道跟随 `soil_heat_method`。

### 11.9 通量与输出

场景通量为各组分按面积或叶面积权重之和：

$$
H_{scene} = \sum_i f_i H_i,\ \lambda E_{scene} = \sum_i f_i (\lambda E)_i,\ G = G_{soil}
$$

输出每步的通量、组分温度、阻抗、稳定度（L、u*、ζ）、闭合标志与残差、各组分的受照面积比（`frac_*_sunlit`）以及诊断量。未闭合的步照常推进热历史，并记录 `closed=False` 与残差大小。

## 12 数值实现

| 项 | 取值 |
|---|---|
| 叶片相函数积分 | 叶倾角 64 点 Gauss（μ 在 [0,1] 均匀） |
| 节点半球积分 | 32 点 Gauss（μ）× 24 方位（行结构关于法向对称，只算一个象限） |
| 地形面元 | 3 个高度带 × 10 个面元方位；锥面半球 9 点 |
| 城市视因子 | 20 点 μ × 72 方位 |
| 冠层节点 | 均质 hom3（土壤 + 2 层）、垄行 row3（覆盖土壤 + 开放土壤 + 1 层）、树冠 crown2（1 个土壤 + 1 层） |
| 光谱波段 | 短波 2001 个 1 nm 波段（0.4–2.4 µm）；方向模拟按输入的波长列表插值 |
| 能量平衡判据 | 最大组分残差 ≤ 5 W m⁻²；阵风一致性残差 ≤ 1e−4 W m⁻² |
| 稳定度判据 | 见 §11.5，外层内部用容差的一半 |
| 稳定度上限 | ζ ≤ 1（强稳定截断），L 的下限由该截断与参考高度决定 |

## 13 一致性检查与已知限制

以下项目由测试集（`tests/`）与物理检查脚本验证。

**已验证的一致性**

- 平地与城市裸地的热红外与短波满足基尔霍夫互易：ε_eff + HDRF(ρ = 1−ε, τ = 0) = 1.00000。
- 垄行闭式间隙与截面数值积分一致（10 组角度与结构，差 < 1e-5）。
- 平地日照土壤的联合间隙对太阳/观测对称（0.0154 / 0.0154）。
- TiRT 与 TIRTEB 的方向观测一致：16 个场景逐方向最大差 2e-16。

**本轮已修复**（2026-10-27）

1. **复合地形裸地的短波反射率曾恒等于土壤反射率**。原因是裸地在 sunlit=False 时仍返回"日照土壤"，地形阴影与面元入射随之丢失。现在裸地遵守光照状态，复合地形裸地按地形比例计算：

$$
R = (1-f_d) r_s \frac{P_{ss}}{P_{ss}+P_{sh}} + f_d r_s
$$

0.86 µm 实测：天底 0.2367，朝太阳 30° 0.2500，背太阳 60° 0.1998（修复前全部为 0.250000）；全散射天空时回到 r_s = 0.25。
2. **短波复合地形的地形参数**现在与热红外一致：缺省 20 m / 10 m / 0.001，支持 `10;20` 分号多层写法，层数按最长列表补齐。
3. **裸地复合地形热红外的直射几何**已改正：锥体截面半径改回 r(1−z/h)（原为 rz/h，随高度增大）；日照份额统一取两条射线的联合可见性 exp[−(Λ_vs−Λ_v)]，取代原先 (1+cos γ)/2 的启发式，且在所有太阳高度下计入他山遮挡。
4. **默认路径上的无效计算已删除**：平地垄行/树冠不再先算一遍旧组分比例再丢弃，通用分支的组分比例只在实际使用处计算。

**已知限制**（属当前闭合方案的固有性质，修改即改变闭合设计）

5. **有植被的复合地形与城市不满足基尔霍夫互易**：ε_eff + HDRF(ρ = 1−ε) 为 1.0026–1.0034（地形）与 1.0035–1.0045（城市），约相当于 0.2–0.3 K。
6. **光照叶片近似导致太阳/观测不互易**：LAI 2、太阳 40° 与观测 70° 对调后，均质 BRF 为 0.378/0.434，垄行 0.390/0.434，树冠 0.319/0.310；来源是 4.5 节"上层可见叶全部按日照计"的假设。
7. **单坡近地平线时 BRF 不随入射余弦线性衰减**：坡 40°、纯直射下 BRF 为 sza 49.9° → 0.0526、49.99° → 0.0053、越过坡面地平线后为 0；BRF/cos i 在最后 0.1° 内从 49 升到 3115。原因是一次散射项不随局地入射余弦消失，且局地太阳天顶角被截断到 89.9°。
8. **树冠热点深度尺度不一致**：土壤用 (4/3)H，日照叶用 2H，同一 hspot 下相关权相差 1.5 倍。
9. **两个波段的坡面背面判据不同**：热红外在相对坡面法线天顶角 > 89° 时省略并返回空值，短波只要求 cos i ≥ 0，同一方向可能短波有值而热红外为空值。
10. **planck 把小于 100 的温度当作摄氏度**（自动加 273.15），输入低于 100 K 会被静默误解释；实际地表温度不会进入该区间。
11. **TIRTEB 的方向反射率与能量平衡反照率不完全自洽**：能量平衡用节点吸收预算（一次逃逸按独立间隙计），方向反射率用四组分一次散射（含热点相关），二者半球积分差 0.016–0.033（近红外、LAI 2）。能量平衡不受影响，但两者不能互相校验。
12. **弱风稳定层结下的交换受经验限幅支配**：min_wind_speed（默认 0.05 m/s）与每段 400 s/m 的阻抗上限在傍晚低风速时起主导作用。默认算例（第 214 天 18:20，风速 0.21 m/s）叶片阻抗被限到 479 s/m，约为正午的 15 倍。

**近似而非缺陷**

- 节点多次散射对直射的首次碰撞也用漫射型后向/前向系数 (2ρ+τ)/3、(ρ+2τ)/3，未使用与太阳天顶角相关的 SAIL 系数。
- 城市墙体受照为整面 0/1 判定；城市街道与地形的两影重叠分别用 max+min·Ω 与 Λ_v+Λ_sΩ 两种闭式。
- 热红外节点散射令叶片 τ = 0、ρ = 1−ε_l，叶片透过率在热红外被并入反射。
