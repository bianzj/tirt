
from old.semiphysical.crown import *
from old.hotspot import *
from matplotlib import pyplot as plt

'''structural variables'''
lai = 5.0
hspot = 0.1
stand = 0.02
radi_horizontal = 1
radi_vertical = 3
'''optical and thermal variables'''
wavelength = 10.5
emissivity_leaf = 0.985
emissivity_soil = 0.955
temperature_leaf_sunlit = 303
temperature_leaf_shaded = 300
temperature_soil_sunlit = 315
temperature_soil_shaded = 305
'''view and solar geometry'''
sza = 30
vza = np.hstack([np.linspace(50,1,50),np.linspace(0,50,51)])
raa = np.hstack([np.repeat(0,51),np.repeat(180,50)])
#raa = 0

'''1.initial'''
crown = Crown()
'''2.set variables'''
crown.set_structure(lai, hspot,stand,radi_horizontal,radi_vertical)
crown.set_optical(wavelength, emissivity_soil, emissivity_leaf)
crown.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
crown.set_angle(vza, sza, raa)
''''3.run'''

Rad = crown.run(ifradiance=1)

vza[raa == 180] = vza[raa == 180]* -1
plt.plot(vza,Rad)
plt.xlabel('VZA ($\circ$)')
plt.ylabel('Brightness Temperature (K)')
plt.show()