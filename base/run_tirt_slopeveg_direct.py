import matplotlib.pyplot as plt
import numpy as np

from base.crown import *

def slope1(tip, phip, tsp, phsp):
    PI = 3.1415926
    # % if (tsp < 0.01)
    #     % sprintf('%s', 'tsp<0.01,returns!');
    # % tip1 = tip;
    # % phip1 = 0;
    # % return;
    # % end

    rd = np.pi/180.0
    sints = np.sin(tsp*rd)
    costs = np.cos(tsp*rd)
    #
    temp = phip - phsp
    cosphi = np.cos(temp*rd)
    sinphi = np.sin(temp*rd)
    #
    sinti = np.sin(tip*rd)
    costi = np.cos(tip*rd)
    #
    x = cosphi * sinti * costs - costi * sints
    y = sinti * sinphi
    z = sinti * cosphi * sints + costs * costi
    r = x * x + y * y
    tip1 = np.arccos(z)
    ind = tip1 > PI/2.0
    tip1[ind] = PI - tip1[ind]
    # if tip1 > PI / 2:
    #     tip1 = PI - tip1


    phip1 = np.zeros_like(r)
    ind  = r!=0
    phip1[ind] = np.arcsin(y[ind] / np.sqrt(r[ind]))

    # if (x < 0):
    #     if (y > 0):
    #         phip1 = 3.14 - phip1
    #     else:
    #         phip1 = -3.14 - phip1
    #
    ind = (x<0)*(y>0)
    phip1[ind] =  3.14 - phip1[ind]
    ind = (x<0)*(y<=0)
    phip1[ind] = -3.14 - phip1[ind]
    ind = r<0.001
    phip1[ind] = 0.0000

    tip1 = tip1/rd
    phip1 = phip1/rd
    return tip1,phip1



def slope0(tip, phip, tsp, phsp):
    PI = 3.1415926
    # % if (tsp < 0.01)
    #     % sprintf('%s', 'tsp<0.01,returns!');
    # % tip1 = tip;
    # % phip1 = 0;
    # % return;
    # % end
    rd = np.pi/180.0
    sints = np.sin(tsp*rd)
    costs = np.cos(tsp*rd)
    #
    temp = phip - phsp
    cosphi = np.cos(temp*rd)
    sinphi = np.sin(temp*rd)
    #
    sinti = np.sin(tip*rd)
    costi = np.cos(tip*rd)
    #
    x = cosphi * sinti * costs - costi * sints
    y = sinti * sinphi
    z = sinti * cosphi * sints + costs * costi
    r = x * x + y * y
    tip1 = np.arccos(z)
    if tip1 > PI / 2:
        tip1 = PI - tip1

    if r <0.0001:
        phip1 = 0
    else:
        phip1 = np.arcsin(y / np.sqrt(r))
        if (x < 0):
            if (y > 0):
                phip1 = 3.14 - phip1
            else:
                phip1 = -3.14 - phip1
        tip1 = tip1/rd
        phip1 = phip1/rd
    return tip1,phip1



def sample_crown():
    '''structural variables'''
    lai = 5.0
    hspot = 0.1
    stand = 0.02
    radi_horizontal = 1
    radi_vertical = 3
    '''optical and thermal variables'''
    wavelength = 10.5
    emissivity_leaf = 0.975
    emissivity_soil = 0.975
    temperature_leaf_sunlit = 22+273.15
    temperature_leaf_shaded = 16+273.15
    temperature_soil_sunlit = 30+273.15
    temperature_soil_shaded = 16+273.15
    '''view and solar geometry'''
    sza = 30
    saa = 0
    vaa = 0
    vza = np.hstack([np.linspace(50, 1, 50), np.linspace(0, 50, 51)])
    raa = np.hstack([np.repeat(0, 51), np.repeat(180, 50)])
    # raa = 0

    '''1.initial'''
    crown = Crown()

    '''2.set variables'''
    crown.set_structure(lai, hspot, stand, radi_horizontal, radi_vertical)
    crown.set_optical(wavelength, emissivity_soil, emissivity_leaf)
    crown.set_thermal(temperature_soil_sunlit,
                      temperature_soil_shaded,
                      temperature_leaf_sunlit,
                      temperature_leaf_shaded)
    crown.set_angle(vza, sza, raa)
    ''''3.run'''

    BT1 = crown.run()


    '''1.initial'''

    pza = 30
    paa = 0
    rd = np.pi/180.0
    '''判断是否是被太阳挡住'''
    '''观测天顶角的等效角度'''
    '''LAI和密度的相应变化'''

    szanew,saanew = slope0(sza,saa,pza,paa)
    vzanew,vaanew = slope1(vza,raa,pza,paa)



    # lainew = lai*cosang
    standnew = stand*np.cos(pza*rd)
    hspotnew = hspot/np.cos(pza*rd)
    # standnew = stand
    # hspotnew = hspot


    crown = Crown()

    '''2.set variables'''
    crown.set_structure(lai, hspotnew, standnew, radi_horizontal, radi_vertical)
    crown.set_optical(wavelength, emissivity_soil, emissivity_leaf)
    crown.set_thermal(temperature_soil_sunlit, temperature_soil_shaded, temperature_leaf_sunlit,
                      temperature_leaf_shaded)
    crown.set_angle(vzanew, szanew, vaanew-saanew)

    ''''3.run'''

    BT2 = crown.run()



    vza[raa == 180] = vza[raa == 180] * -1
    plt.plot(vza, BT1)
    plt.plot(vza, BT2)
    plt.xlabel('VZA ($\\circ$)')
    plt.ylabel('Brightness Temperature (K)')
    plt.title('crown')
    plt.legend(['Veg','Slope'])
    plt.show()
    pass


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    sample_crown()
    # sample_hom()
    # sample_row()


