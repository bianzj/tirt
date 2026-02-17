from old.semiphysical.row import *
from old.physical.row_voxel import *

from old.scatter import *
from old.gap import *


'''structural variables'''
lai = 1.5
hspot = 0.15
row_width = 0.47
row_blank = 0.03
row_height = 1.0
'''optical and thermal variables'''
wavelength = 10.5
emissivity_leaf = 0.985
emissivity_soil = 0.955
temperature_leaf_sunlit = 303
temperature_leaf_shaded = 300
temperature_soil_sunlit = 320
temperature_soil_shaded = 305
'''view and solar geometry'''
vza = np.hstack([np.linspace(50,1,50),np.linspace(0,50,51)])
vaa = np.hstack([np.repeat(0,51),np.repeat(180,50)])
vaa = vaa+0

# vza = 45
# vaa = 0
sza = 25
raa = 0
saa = 0

vsa = np.abs(vaa - saa)

row_voxel = Row_Voxel()
row = Row()


row_voxel.set_structure(lai, hspot, row_width, row_blank, row_height, number_voxel_spacing=50)
row_voxel.set_optical(wavelength, emissivity_soil, emissivity_leaf)
row_voxel.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
row_voxel.set_angle(vza, sza, vaa, saa, raa)

row.set_structure(lai, hspot, row_width, row_blank, row_height)
row.set_optical(wavelength, emissivity_soil, emissivity_leaf)
row.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
row.set_angle(vza, sza, vaa, saa, raa)
''''3.run'''
BT_voxel = row_voxel.run()
BT_hom = row.run()
vza[vsa > 90] = vza[vsa > 90]* -1
plt.plot(BT_voxel,'o-')
plt.plot(BT_hom,'-')
plt.ylim([298,315])
plt.legend(['voxel-based','analytical'])
plt.show()

