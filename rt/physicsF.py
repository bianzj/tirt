
import numpy as np

DEFAULT_WAVELENGTH = 10.5


def planck(wavelength, Ts=None):
    c1 = 11910.439340652
    c2 = 14388.291040407
    if Ts is None:
        Ts = wavelength
        wavelength = DEFAULT_WAVELENGTH

    Ts = np.asarray(Ts, dtype=np.float64)
    wavelength = np.asarray(wavelength, dtype=np.float64)

    if Ts.ndim == 0:
        Ts = Ts.item()
        if (Ts < 100): Ts = Ts + 273.15
        rad = c1 / (np.power(wavelength, 5) * (np.exp(c2 / Ts / wavelength) - 1)) * 10000
    else:
        Ts = Ts.copy()
        Ts[Ts < 100] = Ts[Ts < 100] + 273.15
        rad = c1 / (np.power(wavelength, 5) * (np.exp(c2 / Ts / wavelength) - 1)) * 10000
    return rad




def inv_planck(wavelength, rad=None):
    if rad is None:
        rad = wavelength
        wavelength = DEFAULT_WAVELENGTH
    elif (np.size(wavelength) > 1 and np.size(rad) == 1) or (np.size(wavelength) == 1 and np.size(rad) == 1 and wavelength > 100 and rad <= 100):
        wavelength, rad = rad, wavelength

    c1 = 11910.439340652 * 10000
    c2 = 14388.291040407
    temp = c1 / (rad * np.power((wavelength), 5)) + 1
    Ts = c2 / (wavelength * np.log(temp))
    return Ts




def leaf_inclination_distribution_function(ala):
    '''
    叶倾角分布函数，通过度数ala计算2.5,7.5.。。87.5度的分布概率，总概率的和为1
    :param ala: 角度，度数
    :return: 中间度数，出现概率
    '''

    ### 计算角度步长的中间代表步长

    n = 18
    tx2 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85]
    tx1 = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]
    tx2 = np.asarray(tx2)
    tx1 = np.asarray(tx1)
    x = (tx2 + tx1) / 2.0
    tl1 = tx1 * (np.pi / 180.0)
    tl2 = tx2 * (np.pi / 180.0)
    excent = np.exp(-1.6184e-5 * np.power(ala, 3.) + 2.1145e-3 * np.power(ala, 2.0) - 1.2390e-1 * ala + 3.2491)

    freq = np.zeros(n)
    for i in range(n):
        x1 = excent / (np.sqrt(1.0 + excent * excent * np.tan(tl1[i]) * np.tan(tl1[i])))
        x2 = excent / (np.sqrt(1.0 + excent * excent * np.tan(tl2[i]) * np.tan(tl2[i])))
        if (excent == 1):
            freq[i] = abs(np.cos(tl1[i]) - np.cos(tl2[i]))
        else:
            alpha = excent / np.sqrt(abs(1.0 - excent * excent))
            alpha2 = alpha * alpha
            x12 = x1 * x1
            x22 = x2 * x2
            if (excent > 1):
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
    return x, freq0



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

def ellipsoid_grid_crown_FRT(ncub=45):
    nface = 20
    pi = np.pi
    vol = pi * 4.0 / 3.0

    xtst = np.zeros(ncub)
    ytst = np.zeros(ncub)
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
    for itst in range(2, 14):
        atst[itst - 1] = vol * (491691.0 + 54101.0 * np.sqrt(31.0)) / 21.0924e6
    for itst in range(14, 26):
        atst[itst - 1] = vol * (491691.0 - 54101.0 * np.sqrt(31.0)) / 21.0924e6
    for itst in range(26, 46):
        atst[itst - 1] = vol * 1331.0 / 68.04e3

    alph = np.sqrt(81.0 - 6.0 * np.sqrt(31.0)) / 11.0
    beta = np.sqrt(81.0 + 6.0 * np.sqrt(31.0)) / 11.0
    gamm = 3.0 / np.sqrt(11.0)

    # ceneter of sphere
    xtst[0] = 0.0
    ytst[0] = 0.0
    ztst[0] = 0.0

    # *verteces of icosahedron(12)
    xi[1] = 0.
    xi[2] = 0.
    yi[1] = 0.
    yi[2] = 0.
    zi[1] = 1.
    zi[2] = -1.

    for i in range(4, 9):
        xi[i - 1] = np.cos((i - 4.0) * 2.0 * pi / 5.0) * 2.0 / np.sqrt(5.0)
        xi[i - 1 + 5] = np.cos((2.0 * (i - 4.0) + 1.0) * 2.0 * pi / 5.0) * 2.0 / np.sqrt(5.0)
        yi[i - 1] = np.sin((i - 4.0) * 2.0 * pi / 5.0) * 2.0 / np.sqrt(5.0)
        yi[i - 1 + 5] = np.sin((2.0 * (i - 4.0) + 1.0) * 2.0 * pi / 5.0) * 2.0 / np.sqrt(5.0)
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
    xz = xi[i1[i] - 1] + xi[i2[i] - 1] + xi[i3[i] - 1]
    yz = yi[i1[i] - 1] + yi[i2[i] - 1] + yi[i3[i] - 1]
    zz = zi[i1[i] - 1] + zi[i2[i] - 1] + zi[i3[i] - 1]

    rz = np.sqrt(xz * xz + yz * yz + zz * zz)
    xtst[26 - 1] = 1.0 * xz / rz * gamm
    ytst[26 - 1] = 1.0 * yz / rz * gamm
    ztst[26 - 1] = 1.0 * zz / rz * gamm
    for i in range(1, nface):
        xtst[i + 25] = (xi[i1[i] - 1] + xi[i2[i] - 1] + xi[i3[i] - 1]) / rz * gamm
        ytst[i + 25] = (yi[i1[i] - 1] + yi[i2[i] - 1] + yi[i3[i] - 1]) / rz * gamm
        ztst[i + 25] = (zi[i1[i] - 1] + zi[i2[i] - 1] + zi[i3[i] - 1]) / rz * gamm
    ellipsoid = np.zeros([4, ncub])
    for k in range(ncub):
        ellipsoid[:, k] = [xtst[k], ytst[k], ztst[k], atst[k]]
    return ellipsoid, ncub

def slope1(vza, vaa, pza, paa):
    '''

    :param vza:vza
    :param vaa: vaa
    :param pza: pza
    :param paa: paa
    :return:
    '''
    PI = 3.1415926
    # % if (tsp < 0.01)
    #     % sprintf('%s', 'tsp<0.01,returns!');
    # % tip1 = tip;
    # % phip1 = 0;
    # % return;
    # % end

    rd = np.pi/180.0
    sthetp = np.sin(pza * rd)
    cthetp = np.cos(pza * rd)
    #
    temp = (vaa - paa)

    cosphi = np.cos(temp*rd)
    sinphi = np.sin(temp*rd)
    #
    sthetv = np.sin(vza * rd)
    cthetv = np.cos(vza * rd)
    #
    x = cosphi * sthetv * cthetp - cthetv * sthetp
    y = sthetv * sinphi
    z = sthetv * cosphi * sthetp + cthetp * cthetv
    r = x * x + y * y

    ###--------------------------
    ### cos(theta) = z
    ###--------------------------
    z[z>0.9999] = 0.9999
    tip1 = np.arccos(z)
    ind = tip1 > PI/2.0
    tip1[ind] = PI - tip1[ind]
    # if tip1 > PI / 2:
    #     tip1 = PI - tip1

    phip1 = np.arcsin(y / np.sqrt(r))
    # if (x < 0):
    #     if (y > 0):
    #         phip1 = 3.14 - phip1
    #     else:
    #         phip1 = -3.14 - phip1
    #
    ind = (x<0)*(y>0)
    phip1[ind] =  3.14159 - phip1[ind]
    ind = (x<0)*(y<=0)
    phip1[ind] = -3.14159 - phip1[ind]
    ind = (r)<0.00001
    phip1[ind] = 0.0000

    # ind = phip1 < 0
    # phip1[ind] = phip1[ind] + 3.1415

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



