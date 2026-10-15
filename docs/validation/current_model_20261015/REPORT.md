# 当前TIRT：全部场景光学与热红外结果

按原GitHub测试顺序重新计算10组场景；复合地形替换为已经建立的100 m高、100 m半径山体场景，森林球冠半轴2 m、冠底3 m。原20 m小丘陵和0.9 m叶片森林不参与当前展示。

每个场景依次展示10.5 μm热红外、670 nm红光、860 nm近红外。SZA=30°、SAA=0°；VZA=0–75°、间隔5°；VAA=0/90/180/270°。天底仅存一次，每波段61个方向。每张图的太阳主平面和垂直平面使用相同纵轴范围。

复合山地森林暂定仅统一交换闭合：光学与热红外共用整冠几何和固定4状态交换核心；外热点维持R/hcr，未叠加冠心高度修正。其他场景沿用各自当前模型分支和既有输入。

热红外使用既有固定组分温度与发射率；本次是方向辐射计算，没有重新求解日变化能量平衡。城市森林热红外现已接入整冠内／外双热点，光学原有双热点继续使用；两者默认内宽hspot、外宽R/hcr。建筑遮挡保持城市原有计算，多次散射仍使用城市自身闭合，并未套用山地交换矩阵。

当前TIRT输出来自本次全场景重算。图中蓝线为当前TIRT，橙色圆点曲线为同场景已有HiSTREAM参考；未重新运行HiSTREAM或改变已有几何资产。100 m山地使用对应的大山体参考，其余9场景使用原GitHub场景参考。

光学参考沿用512²、256样本、depth32结果。100 m山地热红外沿用32样本、depth4结果，其中20°/90°使用此前1024²有效重算，其余60方向为512²。原GitHub热红外沿用冻结CSV；参考源和SHA256见reference_manifest.json。

光学沿用既有匹配材料：叶片ρ=0.06/0.45、τ=0.03/0.42，土壤ρ=0.18/0.30，屋顶ρ=0.22/0.35，墙面ρ=0.30/0.45；红光/近红外依次列出。太阳直射照明，天空漫射比例0，内部散射仍计算。单坡继续采用原测试局部坐标和局部入射归一化。

输出1830条波段—方向记录，其中1806条有效预测；单坡背面观测明确标记slope_backface，图上留空。

| 顺序 | 场景 | 当前三波段结果 |
|---|---|---|
| 1 | 平地裸土 | [图](figures/01_T1_plane_bare_all.png) · [输入](inputs/01_T1_plane_bare.json) |
| 2 | 平地均质冠层（LAI 1.5） | [图](figures/02_T2_plane_hom_all.png) · [输入](inputs/02_T2_plane_hom.json) |
| 3 | 平地离散森林（叶片 0.15 m） | [图](figures/04_T3_plane_crown_leaf0.15_all.png) · [输入](inputs/04_T3_plane_crown_leaf0.15.json) |
| 4 | 垄行作物（垄向 0°） | [图](figures/05_T5_row_az0_all.png) · [输入](inputs/05_T5_row_az0.json) |
| 5 | 垄行作物（垄向 90°） | [图](figures/06_T5_row_az90_all.png) · [输入](inputs/06_T5_row_az90.json) |
| 6 | 单坡 30° 朝东 + 均质冠层 | [图](figures/07_S1_slope_hom_all.png) · [输入](inputs/07_S1_slope_hom.json) |
| 7 | 单坡 30° 朝东 + 森林（叶片 0.15 m） | [图](figures/08_S1_slope_crown_all.png) · [输入](inputs/08_S1_slope_crown.json) |
| 8 | 复合山地＋森林（山高100 m、半径100 m） | [图](figures/09_TV_terrain_crown_all.png) · [输入](inputs/09_TV_terrain_crown.json) |
| 9 | 城市（无植被） | [图](figures/10_T4_urban_bare_all.png) · [输入](inputs/10_T4_urban_bare.json) |
| 10 | 城市 + 森林（叶片 0.15 m，31.5 °C） | [图](figures/11_UV_urban_crown_all.png) · [输入](inputs/11_UV_urban_crown.json) |

## 数据与复现

- all_results.csv：全部重新计算的方向结果，热红外value为K，光学value为BRF；热红外图中换算为°C。
- result_ranges.csv：各场景三波段取值范围；timings.csv：本次顺序运行耗时（含缓存复用）。
- histream_reference.csv、reference_comparison.csv、accuracy.csv：同方向参考、配对差值及误差统计；热红外误差单位K，光学为BRF。
- inputs/：实际使用的完整配置；仓库当前源码、manifest.json：本次模型代码副本、源码和输入哈希。
- verification.json：有限值、Planck反算、选定方案复现与坐标检查。

```powershell
python scripts/validate_current_model.py --plots
```
