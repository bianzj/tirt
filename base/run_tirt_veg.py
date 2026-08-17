# This is a sample Python script.

from base.crown import *
from base.row import *
from base.hom import *

def sample_hom():
    '''structural variables'''
    lai = 1.5
    hspot = 0.15
    '''optical and thermal variables'''
    wavelength = 10.5
    emissivity_leaf = 0.985
    emissivity_soil = 0.905
    temperature_leaf_sunlit = 303
    temperature_leaf_shaded = 300
    temperature_soil_sunlit = 320
    temperature_soil_shaded = 305
    '''view and solar geometry'''
    sza = 30
    vza = np.hstack([np.linspace(50, 1, 50), np.linspace(0, 50, 51)])
    raa = np.hstack([np.repeat(0, 51), np.repeat(180, 50)])
    # raa = 0

    '''1.initial;'''
    hom = Hom()
    '''2.set variables;'''
    hom.set_structure(lai, hspot)
    hom.set_optical(wavelength, emissivity_soil, emissivity_leaf)
    hom.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
    hom.set_angle(vza, sza, raa)
    ''''3.run'''
    BT = hom.run()

    vza[raa == 180] = vza[raa == 180] * -1
    plt.plot(vza, BT)
    plt.xlabel('VZA ($\\circ$)')
    plt.ylabel('Brightness Temperature (K)')
    plt.title('hom')
    plt.show()

def sample_row():
    '''structural variables'''
    lai = 1.0
    hspot = 0.5
    row_width = 1
    row_blank = 0.01
    row_height = 1
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
    sza = 30
    raa = 0
    saa = 0

    '''1.initial;'''
    row = Row()
    '''2.set variables;'''
    row.set_structure(lai, hspot, row_width, row_blank, row_height)
    row.set_optical(wavelength, emissivity_soil, emissivity_leaf)
    row.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit, temperature_leaf_shaded)
    row.set_angle(vza, sza, vaa, saa, raa)
    ''''3.run'''
    BT = row.run()

    vza[vaa == 180] = vza[vaa == 180] * -1
    plt.plot(vza, BT)
    plt.xlabel('VZA ($\\circ$)')
    plt.ylabel('Brightness Temperature (K)')
    plt.title('row')
    plt.show()

def sample_crown():
    '''structural variables'''
    lai = 5.0
    hspot = 0.1
    stand = 0.02
    radi_horizontal = 1
    radi_vertical = 3
    '''optical and thermal variables'''
    wavelength = 10.5
    emissivity_leaf = 0.985
    emissivity_soil = 0.955
    temperature_leaf_sunlit = 303
    temperature_leaf_shaded = 300
    temperature_soil_sunlit = 315
    temperature_soil_shaded = 305
    '''view and solar geometry'''
    sza = 30
    vza = np.hstack([np.linspace(50, 1, 50), np.linspace(0, 50, 51)])
    raa = np.hstack([np.repeat(0, 51), np.repeat(180, 50)])
    # raa = 0

    '''1.initial'''
    crown = Crown()
    '''2.set variables'''
    crown.set_structure(lai, hspot, stand, radi_horizontal, radi_vertical)
    crown.set_optical(wavelength, emissivity_soil, emissivity_leaf)
    crown.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit,
                      temperature_leaf_shaded)
    crown.set_angle(vza, sza, raa)
    ''''3.run'''

    BT = crown.run()

    vza[raa == 180] = vza[raa == 180] * -1
    plt.plot(vza, BT)
    plt.xlabel('VZA ($\\circ$)')
    plt.ylabel('Brightness Temperature (K)')
    plt.title('crown')
    plt.show()
    pass


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    sample_crown()
    sample_hom()
    sample_row()


