# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
from base.crown_voxel import *
from base.hom_voxel import *
from base.row_voxel import *

def sample_crown():
    '''structural variables'''
    lai = 3.0

    hspot = 0.10
    stand = 0.07
    hcr = 3.5
    rcr = 2.5

    '''optical and thermal variables'''
    wavelength = 10.5
    emissivity_leaf = 0.975
    emissivity_soil = 0.955
    temperature_leaf_sunlit = 273.15 + 33.5
    temperature_leaf_shaded = 273.15 + 28.9
    temperature_soil_sunlit = 273.15 + 45.8
    temperature_soil_shaded = 273.15 + 30.5
    '''view and solar geometry'''
    sza = 30
    vza = np.hstack([np.linspace(60, 1, 50), np.linspace(0, 60, 51)])
    raa = np.hstack([np.repeat(0, 51), np.repeat(180, 50)])
    # raa = 0

    '''1.initial;'''
    crown_voxel = Crown_Voxel()


    '''2.set variables;'''
    crown_voxel.set_structure(lai, hspot, stand, hcr, rcr, offz=0)
    crown_voxel.set_optical(wavelength, emissivity_soil, emissivity_leaf)
    crown_voxel.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit,temperature_leaf_shaded)
    crown_voxel.set_angle(vza, sza, raa)

    ''''3.run'''
    rad_voxel = crown_voxel.run()
    vza[raa == 180] = vza[raa == 180] * -1
    plt.plot(vza, rad_voxel, 'o-')
    plt.ylim([303, 312])
    plt.legend('voxel-based')
    plt.xlabel('VZA ($\\circ$)')
    plt.ylabel('Brightness Temperature (K)')
    plt.show()

def sample_hom():
    '''structural variables'''
    lai = 3.0
    hspot = 0.35
    '''optical and thermal variables'''
    wavelength = 10.5
    emissivity_leaf = 0.985
    emissivity_soil = 0.955
    temperature_leaf_sunlit = 304
    temperature_leaf_shaded = 300
    temperature_soil_sunlit = 304
    temperature_soil_shaded = 299
    '''view and solar geometry'''
    sza = 30
    vza = np.hstack([np.linspace(50, 1, 50), np.linspace(0, 50, 51)])
    raa = np.hstack([np.repeat(0, 51), np.repeat(180, 50)])
    # raa = 0

    '''1.initial;'''
    hom_voxel = Hom_Voxel()

    '''2.set variables;'''
    hom_voxel.set_structure(lai, hspot, 100)
    hom_voxel.set_optical(wavelength, emissivity_soil, emissivity_leaf)
    hom_voxel.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit,
                          temperature_leaf_shaded)
    hom_voxel.set_angle(vza, sza, raa)
    ''''3.run'''
    BT_voxel = hom_voxel.run()

    vza[raa > 90] = vza[raa > 90] * -1
    plt.plot(vza, BT_voxel, 'o-')
    plt.ylim([301.5, 303.75])
    plt.xlabel('VZA ($\\circ$)')
    plt.ylabel('Brightness Temperature (K)')
    plt.legend('voxel-based')
    plt.show()

def sample_row():
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
    vza = np.hstack([np.linspace(50, 1, 50), np.linspace(0, 50, 51)])
    vaa = np.hstack([np.repeat(0, 51), np.repeat(180, 50)])
    vaa = vaa + 0

    # vza = 45
    # vaa = 0
    sza = 25
    raa = 0 # row azimuth orientation
    saa = 0 # solar azimuth angle

    vsa = np.abs(vaa - saa)

    row_voxel = Row_Voxel()

    row_voxel.set_structure(lai, hspot, row_width, row_blank, row_height, number_voxel_spacing=50)
    row_voxel.set_optical(wavelength, emissivity_soil, emissivity_leaf)
    row_voxel.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit,
                          temperature_leaf_shaded)
    row_voxel.set_angle(vza, sza, vaa, saa, raa)

    ''''3.run'''
    BT_voxel = row_voxel.run()
    vza[vsa > 90] = vza[vsa > 90] * -1
    plt.plot(BT_voxel, 'o-')
    plt.ylim([298, 315])
    plt.legend('voxel-based')
    plt.xlabel('VZA ($\\circ$)')
    plt.ylabel('Brightness Temperature (K)')
    plt.show()


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    sample_crown()
    sample_hom()
    sample_row()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
