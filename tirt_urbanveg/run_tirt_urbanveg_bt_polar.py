# This is a sample Python script.
import matplotlib.pyplot as plt
import numpy as np

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from rt.urban import *
from rt.urban_veg import *
from base.util_file import *
from base.util_plot import *
def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.

'''
被用于高度的分布的模拟
标准化的3连
'''

# Press the green button in the gutter to run the script.
if __name__ == '__main__':


    urban_veg = Urban_Veg()
    urban = Urban()
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


    Estreat = 0.94
    Eroof = 0.92
    Ewall = 0.90
    Eveg = 0.975
    Twall_sunlit = 45 + 273.15
    Twall_shaded = 30+ 273.15
    Tstreet_sunlit = 45+ 273.15
    Tstreet_shaded = 30+ 273.15
    Troof_sunlit = 45+ 273.15
    Troof_shaded = 30+ 273.15
    Tveg = 30 + 273.15

    Bwall_sunlit = planck(wl,Twall_sunlit)
    Bwall_shaded = planck(wl,Twall_shaded)
    Bstreet_sunlit = planck(wl,Tstreet_sunlit)
    Bstreet_shaded = planck(wl,Tstreet_shaded)
    Broof_sunlit = planck(wl,Troof_sunlit)
    Broof_shaded = planck(wl,Troof_shaded)
    Bveg = planck(wl,Tveg)
    B_uv = np.asarray([Bwall_sunlit,Bwall_shaded,Bstreet_sunlit,
                     Bstreet_shaded,Broof_sunlit,Broof_shaded,Bveg])
    B_u = np.asarray([Bwall_sunlit,Bwall_shaded,Bstreet_sunlit,
                     Bstreet_shaded,Broof_sunlit,Broof_shaded])

    ### 计算高度上的差异
    density = 0.003
    # shapes__ = np.asarray([[10, 10, 10, density*0.333,0,90],
    #                      [10, 10, 30, density*0.333,0,90],
    #                      [10, 10, 50, density*0.333,0,90]])
    # shapes_ = np.asarray([[10, 10, 20, density*0.333,0,90],
    #                      [10, 10, 30, density*0.333,0,90],
    #                      [10, 10, 40, density*0.333,0,90]])
    shapes = np.asarray([[10, 10, 30, density,  0,90]])
    forestshape = np.asarray([3.0, 0.01, 3, 3, 0.2])


    urban.set_angular_input(np.abs(vza_),vaa_,sza,saa)
    urban.set_strcutural_input(shapes)
    urban.set_spectral_input(Estreat,Ewall,Eroof)


    emissivity_1 = urban.calculate_effective_component_emissivity(2)
    BB = np.sum(emissivity_1 * B_u,axis=1)
    TB1 = inv_planck(wl,BB)
    # plt.plot(vza_,TB,'o-')
    # plt.xlabel('VZA')
    # plt.ylabel('Brightness Temperatures')

    plt_coutourPolar(vaa_,vza_,TB1,7,61,10)

    # height = 20*0.6+40*0.3+60*0.1
    # shapes = np.asarray([[20,10,30,0.003,0,0+90]])
    urban_veg.set_angular_input(np.abs(vza_),vaa_,sza,saa)
    urban_veg.set_strcutural_input(shapes)
    urban_veg.set_spectral_input(Estreat,Ewall,Eroof,Eveg)
    urban_veg.forestshape = forestshape

    emissivity_1 = urban_veg.calculate_effective_component_emissivity(2)
    BB = np.sum(emissivity_1 * B_uv,axis=1)
    TB2 = inv_planck(wl,BB)
    plt_coutourPolar(vaa_,vza_,TB2-TB1,7,61,10)




# See PyCharm help at https://www.jetbrains.com/help/pycharm/
