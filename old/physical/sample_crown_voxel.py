
from old.physical.crown_voxel import *
from old.semiphysical.crown import *
from old.hotspot import *
from matplotlib import pyplot as plt

'''structural variables'''
lai = 1.5


hspot = 0.1
stand = 0.1
hcr = 2
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
vza = np.hstack([np.linspace(60,1,50),np.linspace(0,60,51)])
raa = np.hstack([np.repeat(0,51),np.repeat(180,50)])
#raa = 0

'''1.initial;'''
crown_voxel = Crown_Voxel()
crown = Crown()

'''2.set variables;'''
crown_voxel.set_structure(lai, hspot, stand, hcr, rcr)
crown_voxel.set_optical(wavelength, emissivity_soil, emissivity_leaf)
crown_voxel.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
crown_voxel.set_angle(vza, sza, raa)
crown.set_structure(lai, hspot,stand,hcr,rcr)
crown.set_optical(wavelength, emissivity_soil, emissivity_leaf)
crown.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
crown.set_angle(vza, sza, raa)
''''3.run'''
rad_voxel =crown_voxel.run(ifradiance=1)
rad = crown.run(ifradiance=1)

vza[raa == 180] = vza[raa == 180]* -1
plt.plot(vza, rad, 'o-')
plt.plot(vza, rad_voxel, 'o-')

plt.legend(['analytical','voxel-based'])
plt.show()

