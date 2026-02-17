
from old.physical.hom_voxel import *
from old.semiphysical.hom import *
from old.hotspot import *
from matplotlib import pyplot as plt

'''structural variables'''
lai = 0.5
hspot = 0.15
'''optical and thermal variables'''
wavelength = 10.5
emissivity_leaf = 0.985
emissivity_soil = 0.955
temperature_leaf_sunlit = 303
temperature_leaf_shaded = 300
temperature_soil_sunlit = 320
temperature_soil_shaded = 305
'''view and solar geometry'''
sza = 25
vza = np.hstack([np.linspace(50,1,50),np.linspace(0,50,51)])
raa = np.hstack([np.repeat(0,51),np.repeat(180,50)])
# raa = 0

'''1.initial;'''
hom_voxel = Hom_Voxel()
hom = Hom()
'''2.set variables;'''
hom_voxel.set_structure(lai, hspot,100)
hom_voxel.set_optical(wavelength, emissivity_soil, emissivity_leaf)
hom_voxel.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
hom_voxel.set_angle(vza, sza, raa)

hom.set_structure(lai,hspot)
hom.set_optical(wavelength,emissivity_soil,emissivity_leaf)
hom.set_thermal(temperature_soil_sunlit,temperature_soil_shaded,temperature_leaf_sunlit,temperature_leaf_shaded)
hom.set_angle(vza,sza,raa)

''''3.run'''
BT_voxel = hom_voxel.run()
BT_hom = hom.run()
vza[raa > 90] = vza[raa > 90]* -1
plt.plot(vza,BT_voxel,'o-')
plt.plot(vza,BT_hom)
plt.ylim([300,315])
plt.xlabel('VZA ($\circ$)')
plt.ylabel('Brightness Temperature (K)')
plt.legend(['voxel-based','analytical'])
plt.show()
