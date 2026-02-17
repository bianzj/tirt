# TiRT 
thermal infrared radiative transfer (TiRT) model


# 典型植被热红外辐射传输建模

>**适用场景**：均质植被、垄行作物和离散森林，植被+城市，植被+地形

>**建模层次**：物理模型、半物理模型

>**建模策略**：体素模型、解析模型

# 改动如下：
1. 去除了核驱动建模部分，完全正向辐射传输模型；
2. 城市代码整合：完成；
   考虑城市建筑异质性模型：tirt_urban;
   考虑城市建筑异质性与街道植被模型：tirt_urbanveg;
3. 山地代码整合：完成；
   考虑单一坡地形与植被模型：tirt_slope;
   考虑复合坡地形与植被模型：tirt_terrain;

地形和城市模型请参考专著：
《复杂地表热红外遥感模型：理论与方法》，作者：卞尊健,肖青,柳钦火

使用说明:

physical,semi-physical 和semi-empirical 目录下分别有均质植被、垄行作物和离散森林的使用案例


0. 调用方法
```
from physical.hom_voxel import *
from semiphysical.hom import *
from hotspot import *
from matplotlib import pyplot as plt
```

1. 设置输入
```
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
vza = np.hstack([np.linspace(50,1,50),np.linspace(0,50,51)])
```
2. 调用类/函数
```
hom_voxel = Hom_Voxel()
hom_voxel.set_structure(lai, hspot,100)
hom_voxel.set_optical(wavelength, emissivity_soil, emissivity_leaf)
hom_voxel.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
hom_voxel.set_angle(vza, sza, raa)
BT_voxel = hom_voxel.run()
```
3. 显示结果
```
plt.plot(BT_voxel)
plt.show()
```


-Email contact: bianzj@aircas.ac.cn

These codes are corresponding to papars as follows and other paper can be found in the reference in these papers:




**Physical**

1. Zunjian Bian, Shengbiao Wu, Jean-Louis Roujean, Biao Cao1, Hua Li, Gaofei Yin, Yongming Du, Qing Xiao, Qinhuo Liu,
A TIR forest reflectance and transmittance (FRT) model for directional temperatures with structural and thermal stratification, Remote Sensing of Environment, 202* 
2. Zunjian Bian, Biao Cao, Hua Li, Yongming Du, Wenjie Fan, Qing Xiao, Qinhuo Liu,
The Effects of Tree Trunks on the Directional Emissivity and Brightness Temperatures of a Leaf-Off Forest Using a Geometric Optical Model, IEEE Transactions on Geoscience and Remote Sensing, 2020a: 1-17

**Semi-physical**

3. Zunjian Bian, Biao Cao, Hua Li, Yongming Du, Jean-Pierre Lagouarde, Qing Xiao,Qinhuo Liu, 
An Analytical Four-component Directional Brightness Temperature Model for Crop and Forest Canopies, Remote Sensing of Environment, 2018, 209(731-746)
4. Zunjian Bian, Qing Xiao, Biao Cao, Yongming Du, Hua Li, Heshun Wang, Qiang Liu,Qinhuo Liu,
Retrieval of Leaf, Sunlit Soil, and Shaded Soil Component Temperatures Using Airborne Thermal Infrared Multiangle Observations, IEEE Transactions on Geoscience and Remote Sensing, 2016, 54(8): 4660-4671

**Semi-empirical**

5. Zunjian Bian, J. L. Roujean, J. P. Lagouarde, Biao Cao, Hua Li, Yongming Du, Qiang Liu, Qing Xiao,Qinhuo Liu, 
A semi-empirical approach for modeling the vegetation thermal infrared directional anisotropy of canopies based on using vegetation indices, ISPRS Journal of Photogrammetry and Remote Sensing, 2020, 160(136-148)
6. Zunjian Bian, Jean-Louis Roujean, Biao Cao, Yongming Du, Hua Li, Philippe Gamet, Junyong Fang, Qing Xiao,Qinhuo Liu,
Modeling the directional anisotropy of fine-scale TIR emissions over tree and crop canopies based on UAV measurements, Remote Sensing of Environment, 2021, 252(112150)
