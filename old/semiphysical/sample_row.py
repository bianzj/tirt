

from old.semiphysical.row import *
from old.hotspot import *
from matplotlib import pyplot as plt

'''structural variables'''
lai = 1.0
hspot = 0.5
row_width = 1
row_spacing = 1.0
row_height = 1
'''optical and thermal variables'''
wavelength = 10.5
emissivity_leaf = 0.985
emissivity_soil = 0.905
temperature_leaf_sunlit = 303
temperature_leaf_shaded = 300
temperature_soil_sunlit = 320
temperature_soil_shaded = 305
'''view and solar geometry'''
vza = np.hstack([np.linspace(50,1,50),np.linspace(0,50,51)])
vaa = np.hstack([np.repeat(0,51),np.repeat(180,50)])
sza = 30
raa = 0
saa = 0

'''1.initial;'''
row = Row()
'''2.set variables;'''
row.set_structure(lai, hspot,row_width,row_spacing,row_height)
row.set_optical(wavelength, emissivity_soil, emissivity_leaf)
row.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
row.set_angle(vza, sza, vaa, saa, raa)
''''3.run'''
BT = row.run()

vza[vaa == 180] = vza[vaa == 180]* -1
plt.plot(vza,BT)
plt.xlabel('VZA ($\circ$)')
plt.ylabel('Brightness Temperature (K)')
plt.show()
