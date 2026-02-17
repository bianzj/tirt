
from old.semiempirical.kdinversion import *
from old.physical.crown_voxel import *
from old.semiempirical.kdmodel import *

### forward simulation

'''structural variables'''
lai = 1.5


hspot = 0.05
stand = 0.1
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
vza = np.hstack([np.linspace(60,1,50),np.linspace(0,60,51)])
vsa = np.hstack([np.repeat(0,51),np.repeat(180,50)])
#raa = 0

'''1.initial;'''
crown_voxel = Crown_Voxel()

'''2.set variables;'''
crown_voxel.set_structure(lai, hspot, stand, hcr, rcr)
crown_voxel.set_optical(wavelength, emissivity_soil, emissivity_leaf)
crown_voxel.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
crown_voxel.set_angle(vza, sza, vsa)

''''3.run'''
bt_voxel =crown_voxel.run()



### fitting
#

coeffs1 = fitting_LSFLiDenseRossThick(bt_voxel,vza,sza,vsa)
bt_sim1 = model_LSFLiDenseRossThick(coeffs1,vza,sza,vsa)

coeffs2 = fitting_VinRL(bt_voxel,vza,sza,vsa)
bt_sim2 = model_VinRL(coeffs2,vza,sza,vsa)


vza[vsa > 90] = vza[vsa > 90]* -1
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

