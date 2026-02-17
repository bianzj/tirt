
from crown_voxel import *
from hotspot import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''structural variables'''
lai = 5.5
hspot = 0.05
stand = 0.08
hcr = 1
rcr = 1
'''optical and thermal variables'''
wavelength = 10.5
emissivity_leaf = 0.985
emissivity_soil = 0.955
temperature_leaf_sunlit = 303
temperature_leaf_shaded = 300
temperature_soil_sunlit = 320
temperature_soil_shaded = 305
'''view and solar geometry'''
sza = 30
vza = [0,70]
raa = [0,0]
#raa = 0

'''1.initial;'''
crown = Crown_Voxel()
'''2.set variables;'''
crown.set_structure(lai, hspot, stand, hcr, rcr)
crown.set_optical(wavelength, emissivity_soil, emissivity_leaf)
crown.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
crown.set_angle(vza, sza, raa)
''''3.run'''
test=crown.run()

plt.plot(test[0],'ko-')
plt.plot(test[1],'ro-')
plt.legend(['voxel-based','analytical'])
plt.show()

