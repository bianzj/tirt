# TIRTEB 理论与方法说明

TIRTEB 是单点、多时相的地表辐射传输与能量平衡模型。每个气象时间步使用太阳几何、地表结构、材料光谱、空气动力学阻抗和植被生理参数，计算组件温度、能量通量及不同方向的观测量。

## 计算流程

~~~text
input_tirteb.csv
  -> input_tirteb.csv / meteo.txt / spectral files / direction files
  -> local solar geometry and surface structure
  -> shortwave and longwave radiation transfer
  -> aerodynamic and surface resistance
  -> photosynthesis or Ball-Berry stomatal resistance
  -> LE, H, G and component energy balance
  -> state update and directional observations
  -> output.csv, observations.csv, intermediate.csv
~~~

相邻时间步继承前一步的土壤温度、组件温度、热历史和深度温度剖面。能量平衡迭代达到收敛条件或 max_iter 上限后进入下一时间步。

## 地表和植被

- plane：平面地表。
- terrain：单坡或复合地形，局地观测角由坡度和坡向转换。
- urban：屋顶、墙面、街道和空气组成的城市场景。
- bare：裸土；hom：均质冠层；row：垄行冠层；crown：树冠/森林冠层。

城市和植被可以组合。城市建筑几何独立决定空气动力学阻抗；植被模型只决定街道冠层的方向组件比例和叶片生理。街道植被覆盖率为 `fveg=urban_vegetation_fraction`，默认 0.30。

## 光谱与观测

可见光和近红外波段使用叶片 PROSPECT 与土壤 BSM 光谱，也支持 measured 文件输入。热红外波段通过叶片反射率、透过率和土壤反射率得到发射率：

~~~text
leaf_emissivity = 1 - leaf_reflectance - leaf_transmittance
soil_emissivity = 1 - soil_reflectance
~~~

热红外组分辐亮度由普朗克函数计算，场景总辐亮度按组件发射率和温度组合得到，再通过反普朗克函数转换为亮温。observations.csv 对小于 2.5 um 的波段保存反射率，对热红外波段默认保存亮温 K。

城市+植被热红外使用两层混合：先用城市几何投影得到墙面、街道和屋顶的方向发射率，再对街道项叠加冠层覆盖。覆盖率为：

~~~text
mu = max(cos(vza), 0.1)
fcover = urban_vegetation_fraction * (1 - exp(-0.5 * LAI / mu))
Lurban+veg = Lurban + fcover * (Pstreet * Lcanopy - Lstreet)
~~~

其中 `Lcanopy` 来自均质、垄行或树冠模型的 turbid 传输辐亮度，`Pstreet` 是城市街道方向可见比例。

观测方向可以使用内置文件、上传文件或直接指定 vza/vaa。观测时间支持全部时间、精确 DOY 和日内小时分数筛选。

## 能量平衡

对土壤、叶片和城市组件，模型计算净辐射、潜热 LE、感热 H 与热储存 G。能量残差可在 intermediate.csv 中检查：

~~~text
residual = net_radiation - LE - H - G
~~~

土壤热量方法支持 direct、surface 和 depth；气孔阻抗支持 Farquhar 和 Ball-Berry；温度修正支持 direct 和 robust。城市建筑组件使用 `urban_rss`，街道冠层叶片使用 Farquhar 或 Ball-Berry 表面阻抗；城市建筑的潜热为 0，植被潜热按覆盖率聚合。

城市空气动力学阻抗采用统一的动量/热量 bulk-transfer 形式。屋顶以气象参考风速和城市冠层上方的位移高度、粗糙度及稳定度计算；墙面和街道使用建筑冠层内的局地街谷风速及由建筑平面尺度确定的街谷粗糙度。街谷阻抗设置相对屋顶的物理下限，以体现建筑遮蔽造成的弱交换。街道内叶片的 `raa_leaf` 在街谷阻抗基础上叠加由 LAI、叶片宽度和局地叶片风速决定的叶片边界层阻抗。

因此城市场景中不再令 `raa_leaf` 直接等于墙面阻抗，且 `raa_soil` 与街道表面使用同一街谷交换阻抗。典型排序为 `raa_roof < raa_street = raa_wall < raa_leaf`；实际数值仍随建筑几何、风速和稳定度变化。

## 输出

- output.csv：每个气象时间步的空气、土壤/叶片温度、LE/H/G、垂直亮温、红光和近红外反射率。
- observations.csv：选定时间、方向和波段的方向性结果。
- intermediate.csv：启用 output_intermediate=yes 后输出的迭代、阻抗、组件比例和能量诊断。

当前实现面向单点时间序列，不是空间栅格模型。论文或方法比较应同时保存有效的 `input_tirteb.csv`、`output.csv`、`observations.csv` 和 `intermediate.csv`。
