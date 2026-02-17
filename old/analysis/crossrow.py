
'''
分析表明垄行和均质场景在温度分布上的差异，当观测结果的数量级或水平相同时候
'''

from old.semiphysical.row import *
from old.semiphysical.hom import *
from old.hotspot import *
from matplotlib import pyplot as plt

'''structural variables'''
lai = 2.5
hspot = 0.15
row_width = 0.5
row_spacing = 0.0
row_height = 2.0
'''optical and thermal variables'''
wavelength = 10.5
emissivity_leaf = 0.985
emissivity_soil = 0.965
temperature_leaf_sunlit = 303
temperature_leaf_shaded = 300
temperature_soil_sunlit = 315
temperature_soil_shaded = 305

dif = 0.00
temperature_leaf_sunlit_ = temperature_leaf_sunlit-dif
temperature_leaf_shaded_ = temperature_leaf_shaded-dif
temperature_soil_sunlit_ = temperature_soil_sunlit-dif
temperature_soil_shaded_ = temperature_soil_shaded-dif
'''view and solar geometry'''

# vza = np.asarray([])
# vaa = np.asarray([])
# for kvaa in range(0,360,30):
#     vza = np.hstack([vza,np.linspace(0,50,11)])
#     vaa = np.hstack([vaa,np.repeat(kvaa,11)])
vaa = [0]
vza = [0]
sza = [20]
raa = np.asarray([])
for kraa in range(0,360,5):
    raa = np.hstack([raa,kraa])
saa = 90

'''1.initial'''
row = Row()
hom = Hom()

'''2.set variables'''
row.set_structure(lai, hspot,row_width,row_spacing,row_height)
row.set_optical(wavelength, emissivity_soil, emissivity_leaf)
row.set_thermal(temperature_soil_sunlit_, temperature_soil_shaded_,
                temperature_leaf_sunlit_, temperature_leaf_shaded_)
row.set_angle(vza, sza, vaa, saa, raa)


hom.set_structure(lai,hspot)
hom.set_optical(wavelength,emissivity_soil,emissivity_leaf)
hom.set_thermal(temperature_soil_sunlit,temperature_soil_shaded,
                temperature_leaf_sunlit,temperature_leaf_shaded)
hom.set_angle(vza,sza,vaa)

# for k in range(90):
#     row.set_angle(vza,sza,vaa,saa,k*1)
#     BT_row = row.run()

'''3.run'''
BT_hom = hom.run()
BT_row = row.run()

plt.plot(BT_hom,'^-')
plt.plot(BT_row,'s-')
plt.show()