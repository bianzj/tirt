
from old.semiempirical.kdinversion import *
from old.physical.row_voxel import *
from old.semiempirical.kdmodel import *

### forward simulation

lai = 1.5
hspot = 0.2
row_width = 0.30
row_blank = 0.20
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
vaa = vaa+90
sza = 30
raa = 0
saa = 90

vsa = np.abs(vaa - saa)


row_voxel = Row_Voxel()
row_voxel.set_structure(lai, hspot, row_width, row_blank, row_height, number_voxel_spacing=50)
row_voxel.set_optical(wavelength, emissivity_soil, emissivity_leaf)
row_voxel.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
row_voxel.set_angle(vza, sza, vaa, saa, raa)
bt_voxel = row_voxel.run()


### fitting
#
coeffs1 = fitting_LSFLiDense(bt_voxel,vza,sza,vsa)
bt_sim1 = model_LSFLiDense(coeffs1,vza,sza,vsa)

# coeffs1 = fitting_LSFLiDenseRossThick(bt_voxel,vza,sza,vsa)
# bt_sim1 = model_LSFLiDenseRossThick(coeffs1,vza,sza,vsa)`

coeffs2 = fitting_VinRL(bt_voxel,vza,sza,vsa)
bt_sim2 = model_VinRL(coeffs2,vza,sza,vsa)



vza[vaa > 90] = vza[vaa > 90]* -1
plt.plot(vza,bt_voxel,'o-')
plt.plot(vza,bt_sim1)
plt.plot(vza,bt_sim2)

dif1 = bt_voxel - bt_sim1
dif2 = bt_voxel - bt_sim2
rmse1 = np.sqrt(np.mean(dif1*dif1))
rmse2 = np.sqrt(np.mean(dif2*dif2))
print(rmse1,rmse2)
plt.ylim([300,315])
plt.xlabel('VZA ($\circ$)')
plt.ylabel('Brightness Temperature (K)')
plt.legend(['voxel-based','$LSFLi_{Dense}$','VinRL'])
plt.show()

