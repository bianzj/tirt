import math

import numpy as np

_EXPINT_NODES = np.arange(1000, 100000) / 1000.0


def _expint(x):
    """C++ Utils::expint numerical integral."""
    if x <= 0.0:
        return float("inf")
    # cumsum keeps the sequential summation order of the C++ loop (np.sum is pairwise).
    return float(np.cumsum(np.exp(-x * _EXPINT_NODES) / _EXPINT_NODES * 0.001)[-1])


def calctav(alpha, nr):
    """C++ SCI::calctav with the same formula and operation order."""
    pi = 3.1415926
    rd = pi / 180.0
    n2 = nr * nr
    np_ = n2 + 1.0
    nm = n2 - 1.0
    a = (nr + 1.0) * (nr + 1.0) / 2.0
    k = -(n2 - 1.0) * (n2 - 1.0) / 4.0
    sa = math.sin(alpha * rd)
    b1 = 0.0
    if alpha != 90.0:
        b1 = math.sqrt((sa * sa - np_ / 2.0) ** 2 + k)
    b2 = sa * sa - np_ / 2.0
    b = b1 - b2
    b3 = b * b * b
    a3 = a * a * a
    ts = (k * k / (6.0 * b3) + k / b - b / 2.0) - (k * k / (6.0 * a3) + k / a - a / 2.0)
    tp1 = -2.0 * n2 * (b - a) / (np_ * np_)
    tp2 = -2.0 * n2 * np_ * math.log(b / a) / (nm * nm)
    tp3 = n2 * (1.0 / b - 1.0 / a) / 2.0
    tp4 = 16.0 * n2 * n2 * (n2 * n2 + 1.0) * math.log((2.0 * np_ * b - nm * nm) / (2.0 * np_ * a - nm * nm)) / (np_ * np_ * np_ * nm * nm)
    tp5 = 16.0 * n2 * n2 * n2 * (1.0 / (2.0 * np_ * b - nm * nm) - 1.0 / (2.0 * np_ * a - nm * nm)) / (np_ * np_ * np_)
    return (ts + tp1 + tp2 + tp3 + tp4 + tp5) / (2.0 * sa * sa)


def _leaf_values(leafbio):
    if isinstance(leafbio, dict):
        return leafbio["Cab"], leafbio["Cw"], leafbio["Cdm"], leafbio["Cs"], leafbio["N"]
    return leafbio[0], leafbio[1], leafbio[2], leafbio[3], leafbio[4]


def fluspect(coeff, leafbio):
    """Calculate leaf reflectance and transmittance as LeafOpt::fluspect."""
    cab, cw, cdm, cs, n = _leaf_values(leafbio)
    reflectance, transmittance = [], []
    for row in coeff[:2001]:
        nr = row[1]
        kab, ks, kw, kdm = row[2], row[4], row[5], row[6]
        kall = (cab * kab + cdm * kdm + cw * kw + cs * ks) / n
        tau = 1.0
        if kall > 0.0:
            tau = (1.0 - kall) * math.exp(-kall) + kall * kall * _expint(kall)

        talf = calctav(59.0, nr)
        ralf = 1.0 - talf
        t12 = calctav(90.0, nr)
        r12 = 1.0 - t12
        t21 = t12 / (nr * nr)
        r21 = 1.0 - t21
        denom = 1.0 - r21 * r21 * tau * tau
        ta = talf * tau * t21 / denom
        ra = ralf + r21 * tau * ta
        tt = t12 * tau * t21 / denom
        r = r12 + r21 * tau * tt
        d = math.sqrt(max((1.0 + r + tt) * (1.0 + r - tt) * (1.0 - r + tt) * (1.0 - r - tt), 0.0))
        rq, tq = r * r, tt * tt
        a = (1.0 + rq - tq + d) / (2.0 * r) if r != 0.0 else 1.0
        b = (1.0 - rq + tq + d) / (2.0 * tt) if tt != 0.0 else 1.0
        b_nm1 = b ** (n - 1.0)
        b_n2_g = b_nm1 * b_nm1
        denom = a * a * b_n2_g - 1.0
        rsub = a * (b_n2_g - 1.0) / denom
        tsub = b_nm1 * (a * a - 1.0) / denom
        if r + tt >= 1.0:
            tsub = tt / (tt + (1.0 - tt) * (n - 1.0))
            rsub = 1.0 - tsub
        denom = 1.0 - rsub * r
        transmittance.append(ta * tsub / denom)
        reflectance.append(ra + ta * rsub * tt / denom)
    return reflectance, transmittance


def leaf_emissivity(reflectance, transmittance):
    return 1.0 - reflectance - transmittance


__all__ = ["calctav", "fluspect", "leaf_emissivity"]
