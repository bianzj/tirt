from base.terrain_vegkd import *
from base.util_plot import *
from base.physicsF import *


def main():
    # t = Terrain()
    # for k in range(2):
    #     d = t.calculate_component_direct_emissivity(k)
    #     s = t.calculate_component_scatter_emissivity(k)
    #     print(d)
    #     print(s)

    terrain = Terrain_Vegkd()
    wl = 10.5
    sza = 30
    saa = 45
    vaa_ = np.asarray([])
    vza_ = np.asarray([])
    vza_temp = np.arange(0,61,10)
    n_temp = np.size(vza_temp)
    for kvaa in range(0,361,10):
        vaa_temp = np.repeat(kvaa,n_temp)
        vza_ = np.hstack([vza_,vza_temp])
        vaa_ = np.hstack([vaa_,vaa_temp])


    Es = 0.94
    Em = 0.92
    Ts_sunlit = 45 + 273.15
    Ts_shaded = 30+ 273.15
    Tm_sunlit = 45 + 273.15
    Tm_shaded = 30+ 273.15


    Bs_sunlit = planck(wl,Ts_sunlit)
    Bs_shaded = planck(wl,Ts_shaded)
    Bm_sunlit = planck(wl,Tm_sunlit)
    Bm_shaded = planck(wl,Tm_shaded)
    B_ = np.asarray([Bs_sunlit,Bs_shaded,Bm_sunlit,
                     Bm_shaded])


    shapes = np.asarray([[10, 3, 0.01]])

    n_shape,n_dim = np.shape(shapes)
    area = 0
    for kshape in range(n_shape):
        area = area + np.pi*shapes[kshape,1]*shapes[kshape,1]*shapes[kshape,2]
    # print('Occupy:',area)

    terrain.set_angular_input(np.abs(vza_),vaa_,sza,saa)
    terrain.set_structural_input(shapes)
    terrain.set_spectral_input(Es,Em)
    emissivity_1 = terrain.calculate_effective_component_emissivity(2)


    BB = np.sum(emissivity_1 * B_,axis=1)
    TB1 = inv_planck(wl,BB)
    # plt.plot(vza_,TB,'o-')
    # plt.xlabel('VZA')
    # plt.ylabel('Brightness Temperatures')

    plt_coutourPolar(vaa_,vza_,TB1-TB1[0],15)


if __name__ == '__main__':
    main()
