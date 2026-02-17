
from old.semiphysical.hom import *
from old.hotspot import *
from matplotlib import pyplot as plt

'''structural variables'''
lai = 1.5
hspot = 0.15
'''optical and thermal variables'''
wavelength = 10.5
emissivity_leaf = 0.985
emissivity_soil = 0.905
temperature_leaf_sunlit = 303
temperature_leaf_shaded = 300
temperature_soil_sunlit = 320
temperature_soil_shaded = 305
'''view and solar geometry'''
sza = 30
vza = np.hstack([np.linspace(50,1,50),np.linspace(0,50,51)])
raa = np.hstack([np.repeat(0,51),np.repeat(180,50)])
# raa = 0

'''1.initial;'''
hom = Hom()
'''2.set variables;'''
hom.set_structure(lai,hspot)
hom.set_optical(wavelength,emissivity_soil,emissivity_leaf)
hom.set_thermal(temperature_soil_sunlit,temperature_soil_shaded,temperature_leaf_sunlit,temperature_leaf_shaded)
hom.set_angle(vza,sza,raa)
''''3.run'''
BT = hom.run()

vza[raa == 180] = vza[raa == 180]* -1
plt.plot(vza,BT)
plt.xlabel('VZA ($\circ$)')
plt.ylabel('Brightness Temperature (K)')
plt.show()
