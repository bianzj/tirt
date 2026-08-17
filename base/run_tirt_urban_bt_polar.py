# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from base.urban import *
from base.util_plot import *
from base.physicsF import *
def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.



# Press the green button in the gutter to run the script.
if __name__ == '__main__':


    urban = Urban()
    wl = 10.5
    sza = 30
    saa = 45
    vaa_ = np.asarray([])
    vza_ = np.asarray([])
    vza_temp = np.arange(0,61,5)
    n_temp = np.size(vza_temp)
    for kvaa in range(0,361,5):
        vaa_temp = np.repeat(kvaa,n_temp)
        vza_ = np.hstack([vza_,vza_temp])
        vaa_ = np.hstack([vaa_,vaa_temp])


    Estreat = 0.94
    Eroof = 0.92
    Ewall = 0.90
    Twall_sunlit = 45 + 273.15
    Twall_shaded = 30+ 273.15
    Tstreet_sunlit = 45+ 273.15
    Tstreet_shaded = 30+ 273.15
    Troof_sunlit = 45+ 273.15
    Troof_shaded = 30+ 273.15

    Bwall_sunlit = planck(wl,Twall_sunlit)
    Bwall_shaded = planck(wl,Twall_shaded)
    Bstreet_sunlit = planck(wl,Tstreet_sunlit)
    Bstreet_shaded = planck(wl,Tstreet_shaded)
    Broof_sunlit = planck(wl,Troof_sunlit)
    Broof_shaded = planck(wl,Troof_shaded)
    B_ = np.asarray([Bwall_sunlit,Bwall_shaded,Bstreet_sunlit,
                     Bstreet_shaded,Broof_sunlit,Broof_shaded])

    # shapes = np.asarray([[10, 10, 10, 0.003*0.333,0,90],
    #                      [10, 10, 30, 0.003*0.333,0,90],
    #                      [10, 10, 50, 0.003*0.333,0,90]])
    # shapes_ = np.asarray([[10, 10, 30, 0.003,  0,90]])

    shapes = np.asarray([[10, 5, 30, 0.003,0,90]])
    shapes_ = np.asarray([[5, 10, 30, 0.003, 0,90]])


    # shapes = np.asarray([[5, 20, 30, 0.003*0.333,0,90],
    #                      [10, 10, 30, 0.003*0.333,0,90],
    #                      [20, 5, 30, 0.003*0.333,0,90]])
    # shapes_ = np.asarray([[10, 10, 30, 0.003,  0,90]])

    # shapes = np.asarray([[10, 10, 30, 0.003*0.333,135+40-90,135+40],
    #                      [10, 10, 30, 0.003*0.333,135-90,135],
    #                      [10, 10, 30, 0.003*0.333,135-40-90,135-40]])
    # shapes_ = np.asarray([[10, 10, 30, 0.003,  135 - 90,135]])




    urban.set_angular_input(np.abs(vza_),vaa_,sza,saa)
    urban.set_structural_input(shapes)
    urban.set_spectral_input(Estreat,Ewall,Eroof)
    # emissivity_1 = tirt_urban.calculate_effective_component_emissivity(1)
    #

    # plt.plot(vza_,emissivity_1,'o-')
    # # plt.plot(vza_,emissivity_2,'k')
    # plt.xlabel('VZA')
    # plt.ylabel('Effective Emissivity')
    # plt.legend(['Wall','Street','Roof'])
    # plt.show()
    # emissivity_11 = np.sum(emissivity_1,axis=1)
    # plt.figure(figsize=[4.5,4])
    # plt.plot(vza_, emissivity_11, 'o-')
    # plt.xlabel('VZA')
    # plt.ylabel('Effective Emissivity')
    # plt.show()
    # print(emissivity_1)


    emissivity_1 = urban.calculate_effective_component_emissivity(2)
    BB = np.sum(emissivity_1 * B_,axis=1)
    TB1 = inv_planck(wl,BB)
    # plt.plot(vza_,TB,'o-')
    # plt.xlabel('VZA')
    # plt.ylabel('Brightness Temperatures')

    plt_coutourPolar(vaa_,vza_,TB1,8)


    # height = 20*0.6+40*0.3+60*0.1
    # shapes = np.asarray([[20,10,30,0.003,0,0+90]])
    urban.set_angular_input(np.abs(vza_),vaa_,sza,saa)
    urban.set_structural_input(shapes_)
    urban.set_spectral_input(Estreat,Ewall,Eroof)

    emissivity_1 = urban.calculate_effective_component_emissivity(2)
    BB = np.sum(emissivity_1 * B_,axis=1)
    TB2 = inv_planck(wl,BB)
    plt_coutourPolar(vaa_,vza_,TB2,8)

    # plt_coutourPolar(vaa_, vza_, TB2-TB1, 8)
    # plt.plot(vza_,TB,'o-')
    # plt.xlabel('VZA')
    # plt.ylabel('Brightness Temperatures')
    # plt.legend(['het','hom'])
    # plt.show()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
