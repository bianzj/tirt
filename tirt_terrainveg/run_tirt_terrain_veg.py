from rt.terrain_veg import *
from base.util_plot import *


def main():
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
    shapes = np.asarray([[10, 10, 0.001],[10, 10, 0.001]])


    forestshape = np.asarray([lai,stand,rcr,hcr,hspot])
    print('terrain:',100*0.002*np.pi)
    print('tree:',np.pi*rcr*rcr*stand)
    terrain = Terrain_Veg()
    crown = Crown()

    vaa_ = np.asarray([])
    vza_ = np.asarray([])
    vza_temp = np.arange(0,66,5)
    n_temp = np.size(vza_temp)
    for kvaa in range(0,361,5):
        vaa_temp = np.repeat(kvaa,n_temp)
        vza_ = np.hstack([vza_,vza_temp])
        vaa_ = np.hstack([vaa_,vaa_temp])


    terrain.forestshape = forestshape
    Bs_sunlit = planck(wl,Ts_sunlit)
    Bs_shaded = planck(wl,Ts_shaded)
    Bv_sunlit = planck(wl,Tv_sunlit)
    Bv_shaded = planck(wl,Tv_shaded)
    B_ = np.asarray([Bs_sunlit,Bs_shaded,Bv_sunlit,Bv_shaded])

    n_shape,n_dim = np.shape(shapes)
    area = 0
    for kshape in range(n_shape):
        area = area + np.pi*shapes[kshape,1]*shapes[kshape,1]*shapes[kshape,2]
    print('Occupy:',area)
    terrain.set_angular_input(np.abs(vza_),vaa_,sza,saa)
    terrain.set_structural_input(shapes)
    terrain.set_spectral_input(Es,El)
    emissivity_1 = terrain.calculate_effective_component_emissivity(2)
    BB = np.sum(emissivity_1 * B_,axis=1)
    # print(BB)
    BTT = inv_planck(wl,BB)
    # plt.plot(vza_,TB,'o-')
    # plt.xlabel('VZA')
    # plt.ylabel('Brightness Temperatures')
    DBTT = BTT - BTT[0]
    plt_coutourPolar(vaa_,vza_,DBTT,7,66)

    '''1.initial'''
    crown = Crown()
    '''2.set variables'''
    crown.set_structure(lai, hspot, stand, rcr, hcr)
    crown.set_optical(wl, Es, El)
    crown.set_thermal(Ts_sunlit, Ts_shaded, Tv_sunlit,Tv_shaded)
    crown.set_angle(vza_, sza, vaa_)
    ''''3.run'''

    BTP = crown.run()
    DBTP = BTP - BTP[0]

    plt_coutourPolar(vaa_,vza_,DBTP,7,66)


    DT = DBTT - DBTP
    plt_coutourPolar(vaa_,vza_,DT,7,66)







if __name__ == '__main__':
    main()
