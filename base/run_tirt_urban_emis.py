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

    sza = 20
    saa = 180
    vza_ = np.arange(0,60,5)
    n_angle = np.size(vza_)
    vaa_ = np.repeat(0,n_angle)
    Estreat = 0.94
    Eroof = 0.92
    Ewall = 0.90

    Twall_sunlit = 45
    Twall_shaded = 30
    Tstreet_sunlit = 45
    Tstreet_shaded = 30
    Troof_sunlit = 45
    Troof_shaded = 30

    Bwall_sunlit = planck(wl,Twall_sunlit)
    Bwall_shaded = planck(wl,Twall_shaded)
    Bstreet_sunlit = planck(wl,Tstreet_sunlit)
    Bstreet_shaded = planck(wl,Tstreet_shaded)
    Broof_sunlit = planck(wl,Troof_sunlit)
    Broof_shaded = planck(wl,Troof_shaded)


    shapes = np.asarray([[10, 10, 10, 0.003*0.5,0,90],
                         [10, 10, 50, 0.003*0.5,0,90]])

    # shapes = np.asarray([[10,10,30,0.003]])

    urban.set_angular_input(vza_,vaa_,sza,saa)
    urban.set_structural_input(shapes)
    urban.set_spectral_input(Estreat,Ewall,Eroof)
    emissivity_1 = urban.calculate_effective_emissivity()

    # height = 20*0.6+40*0.3+60*0.1
    # shapes = np.asarray([[14,7,height,0.003]])

    shapes = np.asarray([[10,10,30,0.003,0,90]])
    urban.set_angular_input(vza_,vaa_,sza,saa)
    urban.set_structural_input(shapes)
    urban.set_spectral_input(Estreat,Ewall,Eroof)
    emissivity_2 = urban.calculate_effective_emissivity()

    plt.plot(vza_,emissivity_1,'r')
    plt.plot(vza_,emissivity_2,'k')
    plt.show()


    a =urban.calculate_direct_emissivity(0,1)
    print(a)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
