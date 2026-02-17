import numpy as np

def planck(Ts,wavelength=10.5):
    c1 = 11910.439340652
    c2 = 14388.291040407

    if isinstance(Ts * 1.0, float):
        if (Ts < 100): Ts = Ts + 273.15
        wavelength = np.float(wavelength)
        Ts = np.float(Ts)
        rad = c1 / (np.power(wavelength, 5) * (np.exp(c2 / Ts / wavelength) - 1)) * 10000
    else:
        Ts[Ts < 100] = Ts[Ts < 100] + 273.15
        wavelength = np.float(wavelength)
        rad = c1 / (np.power(wavelength, 5) * (np.exp(c2 / Ts / wavelength) - 1)) * 10000
    return rad

def inv_planck(rad, wavelength= 10.5):
    c1 = 11910.439340652 * 10000
    c2 = 14388.291040407

    temp = c1 / (rad * np.power((wavelength), 5)) + 1
    Ts = c2 / (wavelength * np.log(temp))
    return Ts

def gap_average(lai):
    m = np.exp(-lai * 0.825)
    return m

def ellipsoid_grid_display(r, c, ng, xg):
    # *****************************************************************************80
    #
    ## ELLIPSOID_GRID_DISPLAY displays grid points inside a ellipsoid.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    11 April 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Parameters:
    #
    #    Input, real R[3], the half axis lengths.
    #
    #    Input, real C[3], the center of the ellipsoid.
    #
    #    Input, integer NG, the number of grid points inside the ellipsoid.
    #
    #    Input, real XYZ[NG,3], the grid point coordinates.
    #
    #    Input, string FILENAME, the name of the plotfile to be created.
    #
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xg[:, 0], xg[:, 1], xg[:, 2], 'b');

    ax.set_xlabel('<---X--->')
    ax.set_ylabel('<---Y--->')
    ax.set_zlabel('<---Z--->')
    ax.set_title('Grid points in ellipsoid')
    ax.grid(True)
    # ax.axis ( 'equal' )
    # plt.savefig(filename)
    plt.show(block=False)
    plt.clf()

    return

def ellipsoid_grid_crown_FRT(ncub = 45):
    nface = 20
    pi = np.pi
    vol = pi * 4.0 / 3.0

    xtst = np.zeros(ncub)
    ytst =  np.zeros(ncub)
    ztst = np.zeros(ncub)
    atst = np.zeros(ncub)
    xi = np.zeros(ncub)
    yi = np.zeros(ncub)
    zi = np.zeros(ncub)
    i1 = np.asarray([2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
          4, 4, 4, 5, 5, 6, 6, 7, 7, 8])
    i2 = np.asarray([4, 4, 5, 6, 7, 9, 9, 10, 11, 12,
          5, 8, 9, 6, 9, 7, 10, 8, 11, 12])
    i3 = np.asarray([5, 8, 6, 7, 8, 10, 13, 11, 12, 13,
          9, 13, 13, 10, 10, 11, 11, 12, 12, 13])
    # i1 = np.asarray(i1) - 1
    # i2 = np.asarray(i2) - 1
    # i3 = np.asarray(i3) - 1
    atst[0] = vol * 2096.0 / 42525.0
    for itst in range(2,14):
        atst[itst - 1] = vol * (491691.0 + 54101.0 * np.sqrt(31.0)) / 21.0924e6
    for itst in range(14,26):
        atst[itst - 1] = vol * (491691.0 - 54101.0 * np.sqrt(31.0)) / 21.0924e6
    for itst in range (26,46):
        atst[itst - 1] = vol * 1331.0 / 68.04e3


    alph = np.sqrt(81.0 - 6.0 * np.sqrt(31.0)) / 11.0
    beta = np.sqrt(81.0 + 6.0 * np.sqrt(31.0)) / 11.0
    gamm = 3.0 / np.sqrt(11.0)

    #ceneter of sphere
    xtst[0] = 0.0
    ytst[0] = 0.0
    ztst[0] = 0.0

    #*verteces of icosahedron(12)
    xi[1] = 0.
    xi[2] = 0.
    yi[1] = 0.
    yi[2] = 0.
    zi[1] = 1.
    zi[2] = -1.


    for i in range(4,9):
        xi[i - 1] = np.cos((i - 4.0) * 2.0 * pi / 5.0) * 2.0 / np.sqrt(5.0)
        xi[i - 1 + 5] = np.cos((2.0 * (i - 4.0) + 1.0) * 2.0 * pi / 5.0)*2.0 / np.sqrt(5.0)
        yi[i - 1] = np.sin((i - 4.0) * 2.0 * pi / 5.0) * 2.0 / np.sqrt(5.0)
        yi[i - 1 + 5] = np.sin((2.0 * (i - 4.0) + 1.0) * 2.0 * pi / 5.0)*2.0 / np.sqrt(5.0)
        zi[i - 1] = 1.0 / np.sqrt(5.0)
        zi[i - 1 + 5] = -1.0 / np.sqrt(5.0)

    for i in range(1, 13):
        xtst[i] = xi[i] * alph
        ytst[i] = yi[i] * alph
        ztst[i] = zi[i] * alph
        xtst[i + 12] = xi[i] * beta
        ytst[i + 12] = yi[i] * beta
        ztst[i + 12] = zi[i] * beta

    # projections of facets centers

    i = 0
    xz = xi[i1[i]-1] + xi[i2[i]-1] + xi[i3[i]-1]
    yz = yi[i1[i]-1] + yi[i2[i]-1] + yi[i3[i]-1]
    zz = zi[i1[i]-1] + zi[i2[i]-1] + zi[i3[i]-1]

    rz = np.sqrt(xz*xz + yz*yz + zz*zz)
    xtst[26 - 1] = 1.0 * xz / rz * gamm
    ytst[26 - 1] = 1.0 * yz / rz * gamm
    ztst[26 - 1] = 1.0 * zz / rz * gamm
    for i in range(1, nface):
        xtst[i + 25] = (xi[i1[i]-1] + xi[i2[i]-1] + xi[i3[i]-1]) / rz * gamm
        ytst[i + 25] = (yi[i1[i]-1] + yi[i2[i]-1] + yi[i3[i]-1]) / rz * gamm
        ztst[i + 25] = (zi[i1[i]-1] + zi[i2[i]-1] + zi[i3[i]-1]) / rz * gamm
    ellipsoid = np.zeros([4,ncub])
    for k in range(ncub):
        ellipsoid[:,k] = [xtst[k],ytst[k],ztst[k],atst[k]]
    return ellipsoid,ncub

def projection_area_stem(z1, z2, dbh, hc, check=0):
    if (check == 0):
        dz = z2 - z1
        dz[dz < 0.001] = 0.001
        return dbh * dz
    ati = [118.9810, -277.5780, 1140.5250, -3037.4870, 4419.6820, -3361.780, 997.6570]
    ht0 = 26.
    dt0 = 30.
    pty = 0.007
    qty = -0.007

    xz0 = 1.30 / hc
    xz1 = z1 / hc
    xz2 = z2 / hc
    eet = pty * (hc - ht0) + qty * (dbh - dt0)

    sum1 = 0.0
    for j in range(1, 7):
        sum1 = sum1 + ati[j - 1] * (np.power(xz0, j - 1))

    f13 = sum1 * (1.0 + eet * (xz0 * xz0 - 0.010))

    sum1 = 0.
    sum2 = 0.
    sum3 = 0.
    sum4 = 0.
    for j in range(1, 7):
        sum1 = sum1 + ati[j - 1] / (j) * (np.power(xz1, j))
        sum2 = sum2 + ati[j - 1] / (j) * (np.power(xz2, j))
        sum3 = sum3 + ati[j - 1] / (j + 3) * np.power(xz1, (j + 3))
        sum4 = sum4 + ati[j - 1] / (j + 3) * np.power(xz2, (j + 3))

    stemArea = (sum2 - sum1) * (1.0 - 0.010 * eet)
    stemArea = stemArea + eet * (sum4 - sum3)
    stemArea = stemArea * hc * dbh / f13
    return stemArea

#######################################################

def leaf_inclination_distribution_function(ala):
    tx2 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85]
    tx1 = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]
    tx2 = np.asarray(tx2)
    tx1 = np.asarray(tx1)
    n = 18
    x = (tx2 + tx1) / 2.0
    tl1 = tx1 * (np.pi / 180.0)
    tl2 = tx2 * (np.pi / 180.0)
    excent = np.exp(-1.6184e-5 * np.power(ala,3.) + 2.1145e-3 * np.power(ala,2.0) - 1.2390e-1 * ala + 3.2491)
    sum0 = 0

    freq = np.zeros(18)
    for i in range(n):
        x1 = excent / (np.sqrt(1.0 + excent*excent * np.tan(tl1[i])*np.tan(tl1[i])))
        x2 = excent / (np.sqrt(1.0 + excent*excent * np.tan(tl2[i])*np.tan(tl2[i])))
        if (excent == 1):
            freq[i] = abs(np.cos(tl1[i]) - np.cos(tl2[i]))
        else:
            alpha = excent / np.sqrt(abs(1.0 - excent*excent))
            alpha2 = alpha * alpha
            x12 = x1*x1
            x22 = x2*x2
            if (excent> 1):
                alpx1 = np.sqrt(alpha2 + x12)
                alpx2 = np.sqrt(alpha2 + x22)
                dum = x1 * alpx1 + alpha2 * np.log(x1 + alpx1)
                freq[i] = abs(dum - (x2 * alpx2 + alpha2 * np.log(x2 + alpx2)))
            else:
                almx1 = np.sqrt(alpha2 - x12)
                almx2 = np.sqrt(alpha2 - x22)
                dum = x1 * almx1 + alpha2 * np.arcsin(x1 / alpha)
                freq[i] = abs(dum - (x2 * almx2 + alpha2 * np.arcsin(x2 / alpha)))


    sum0 = np.sum(freq)
    freq0 = freq / sum0
    return freq0

def atmospheric_correction_leaf(Ts_leaf,Ts_sky,emis_leaf,wl = 10.5):
    if (Ts_leaf < 100):Ts = Ts_leaf+273.15
    if (Ts_sky<100): Ts_sky = Ts_sky+273.15
    rad = planck(Ts_leaf,wl)
    rad_sky = planck(Ts_sky,wl)
    rads = rad-rad_sky*(1-emis_leaf)
    radreal = rads/emis_leaf
    Tsreal = inv_planck(radreal)
    return Tsreal

def atmospheric_correction_soil(Ts_soil,Ts_sky,Ts_leaf,lai,emis_leaf,emis_soil,wl = 10.5):
    if (Ts_soil < 100):Ts_soil = Ts_soil+273.15
    if (Ts_leaf < 100): Ts_leaf = Ts_leaf + 273.15
    if (Ts_sky<100): Ts_sky = Ts_sky+273.15
    rad_leaf = planck(Ts_leaf,wl)
    rad_soil = planck(Ts_soil, wl)
    rad_sky = planck(Ts_sky,wl)
    M = gap_average(lai)
    rad_down = rad_sky*M+rad_leaf*(1-M)*emis_leaf
    rads = rad_soil-rad_down*(1-emis_soil)
    radreal = rads/emis_soil
    Tsreal = inv_planck(radreal)
    return Tsreal

def atmospheric_effect_leaf(Ts_leaf,Ts_sky,emis_leaf,wl = 10.5):
    if (Ts_leaf < 100):Ts = Ts_leaf+273.15
    if (Ts_sky<100): Ts_sky = Ts_sky+273.15
    rad = planck(Ts_leaf,wl)
    rad_sky = planck(Ts_sky,wl)
    rads = rad*emis_leaf+rad_sky*(1-emis_leaf)
    Tsreal = inv_planck(rads)
    return Tsreal

def atmospheric_effect_soil(Ts_soil,Ts_sky,Ts_leaf,lai,emis_leaf,emis_soil,wl = 10.5):
    if (Ts_soil < 100):Ts_soil = Ts_soil+273.15
    if (Ts_leaf < 100): Ts_leaf = Ts_leaf + 273.15
    if (Ts_sky<100): Ts_sky = Ts_sky+273.15
    rad_leaf = planck(Ts_leaf,wl)
    rad_soil = planck(Ts_soil, wl)
    rad_sky = planck(Ts_sky,wl)
    M = gap_average(lai)
    rad_down = rad_sky*M+rad_leaf*(1-M)*emis_leaf
    rads = rad_soil*emis_soil+rad_down*(1-emis_soil)
    Tsreal = inv_planck(rads,wl)
    return Tsreal


