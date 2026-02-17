
from rt.slope_veg import *
from rt.crown import *
from base.util_plot import *

# t = Terrain()
# for k in range(2):
#     d = t.calculate_component_direct_emissivity(k)
#     s = t.calculate_component_scatter_emissivity(k)
#     print(d)
#     print(s)
lai = 3.0
stand = 0.02
hspot = 0.2
rcr = 3
hcr = 3

print(np.pi*rcr*rcr*stand)

Es = 0.955
El = 0.975
Ts_sunlit = 45 + 273.15
Ts_shaded = 30 + 273.15
Tv_sunlit = 33 + 273.15
Tv_shaded = 30 + 273.15
wl = 10.5
sza = 46.1
saa = 0.1
pza = 3.0
paa = 0.0

forestshape = np.asarray([lai,stand,rcr,hcr,hspot])
print('terrain:',100*0.002*np.pi)
print('tree:',np.pi*rcr*rcr*stand)
slope = Slope_Veg()
crown = Crown()

vaa_ = np.asarray([])
vza_ = np.asarray([])
vza_temp = np.arange(0,66,5)
n_temp = np.size(vza_temp)
for kvaa in range(0,361,5):
    vaa_temp = np.repeat(kvaa,n_temp)
    vza_ = np.hstack([vza_,vza_temp])
    vaa_ = np.hstack([vaa_,vaa_temp])


slope.forestshape = forestshape
Bs_sunlit = planck(wl,Ts_sunlit)
Bs_shaded = planck(wl,Ts_shaded)
Bv_sunlit = planck(wl,Tv_sunlit)
Bv_shaded = planck(wl,Tv_shaded)
B_ = np.asarray([Bs_sunlit,Bs_shaded,Bv_sunlit,Bv_shaded])


slope.set_angle(np.abs(vza_), sza, vaa_, saa)
slope.set_optical(wl,Es, El)
slope.set_structure(lai,hspot,stand,hcr,rcr)
slope.set_slope(pza,paa)
slope.set_thermal(Ts_sunlit,Ts_shaded,Tv_sunlit,Tv_shaded)
bt_slopeveg = slope.run()

plt_coutourPolar(vaa_,vza_,bt_slopeveg,7,66)



# #
# # '''1.initial'''
crown = Crown()
'''2.set variables'''
crown.set_structure(lai, hspot, stand, rcr, hcr)
crown.set_optical(wl, Es, El)
crown.set_thermal(Ts_sunlit, Ts_shaded, Tv_sunlit,Tv_shaded)
crown.set_angle(vza_, sza, vaa_)
''''3.run'''

bt_veg = crown.run()
plt_coutourPolar(vaa_,vza_,bt_veg,7,66)
plt.show()

DT = bt_slopeveg - bt_veg
plt_coutourPolar(vaa_,vza_,DT,7,66)
