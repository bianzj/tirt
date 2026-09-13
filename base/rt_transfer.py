import math

import numpy as np

from .gap import gap_probability_hom_analytical
from .hotspot import hotspot_analytical, hotspot_layer, hotspot_vegetation_volume


PI = math.pi
RD = PI / 180.0
SIGMA = 5.67e-8


def jfunc1(k, ell, thickness):
    delta = (k - ell) * thickness
    if delta > 1e-3:
        return (math.exp(-ell * thickness) - math.exp(-k * thickness)) / max(k - ell, 1e-12)
    return 0.5 * thickness * (math.exp(-k * thickness) + math.exp(-ell * thickness)) * (1.0 - delta * delta / 12.0)


def jfunc2(k, ell, thickness):
    return (1.0 - math.exp(-(k + ell) * thickness)) / max(k + ell, 1e-12)


def volscatt(sza, vza, raa, ttl):
    costs, costo = math.cos(RD * sza), math.cos(RD * vza)
    sints, sinto = math.sin(RD * sza), math.sin(RD * vza)
    costl, sintl = math.cos(RD * ttl), math.sin(RD * ttl)
    cosraa = math.cos(RD * raa)
    cs, co = costl * costs, costl * costo
    ss, so = sintl * sints, sintl * sinto
    cosbts = -cs / ss if abs(ss) > 1e-6 else 5.0
    cosbto = -co / so if abs(so) > 1e-6 else 5.0
    if abs(cosbts) < 1.0:
        bts, ds = math.acos(cosbts), ss
    else:
        bts, ds = PI, cs
    chi_s = 2.0 / PI * ((bts - PI / 2.0) * cs + math.sin(bts) * ss)
    if abs(cosbto) < 1.0:
        bto, doo = math.acos(cosbto), so
    elif vza < 90.0:
        bto, doo = PI, co
    else:
        bto, doo = 0.0, -co
    chi_o = 2.0 / PI * ((bto - PI / 2.0) * co + math.sin(bto) * so)
    btran1 = abs(bts - bto)
    btran2 = PI - abs(bts + bto - PI)
    raar = RD * abs(raa)
    if raar <= btran1:
        bt1, bt2, bt3 = raar, btran1, btran2
    elif raar <= btran2:
        bt1, bt2, bt3 = btran1, raar, btran2
    else:
        bt1, bt2, bt3 = btran1, btran2, raar
    t1 = 2.0 * cs * co + ss * so * cosraa
    t2 = math.sin(bt2) * (2.0 * ds * doo + ss * so * math.cos(bt1) * math.cos(bt3)) if bt2 > 0 else 0.0
    frho = max(((PI - bt2) * t1 + t2) / (2.0 * PI * PI), 0.0)
    ftau = max((-bt2 * t1 + t2) / (2.0 * PI * PI), 0.0)
    return chi_s, chi_o, frho, ftau


def cal_lidf(ala):
    tx2 = np.arange(0.0, 90.0, 5.0)
    tx1 = tx2 + 5.0
    excent = math.exp(-1.6184e-5 * ala**3 + 2.1145e-3 * ala**2 - 1.2390e-1 * ala + 3.2491)
    freq = np.zeros(18)
    for i in range(18):
        tl1, tl2 = tx1[i] * RD, tx2[i] * RD
        x1 = excent / math.sqrt(1.0 + excent * excent * math.tan(tl1) ** 2)
        x2 = excent / math.sqrt(1.0 + excent * excent * math.tan(tl2) ** 2)
        if abs(excent - 1.0) < 1e-8:
            freq[i] = abs(math.cos(tl1) - math.cos(tl2))
        else:
            alpha = excent / math.sqrt(abs(1.0 - excent * excent))
            alpha2 = alpha * alpha
            if excent > 1.0:
                a1, a2 = math.sqrt(alpha2 + x1 * x1), math.sqrt(alpha2 + x2 * x2)
                f1 = x1 * a1 + alpha2 * math.log(max(x1 + a1, 1e-12))
                f2 = x2 * a2 + alpha2 * math.log(max(x2 + a2, 1e-12))
            else:
                a1, a2 = math.sqrt(max(alpha2 - x1 * x1, 0.0)), math.sqrt(max(alpha2 - x2 * x2, 0.0))
                f1 = x1 * a1 + alpha2 * math.asin(max(-1.0, min(1.0, x1 / alpha)))
                f2 = x2 * a2 + alpha2 * math.asin(max(-1.0, min(1.0, x2 / alpha)))
            freq[i] = abs(f1 - f2)
    return freq / max(float(freq.sum()), 1e-12)


def _scattering_coefficients(sza, vza, raa, ala):
    cts = max(math.cos(RD * sza), 1e-6)
    cto = max(math.cos(RD * vza), 1e-6)
    ctscto = max(cts * cto, 1e-6)
    freq = cal_lidf(ala)
    ks = ko = bf = sob = sof = 0.0
    for k in range(18):
        ttl = 2.5 + 5.0 * k
        chi_s, chi_o, frho, ftau = volscatt(sza, vza, raa, ttl)
        ks += chi_s / cts * freq[k]
        ko += chi_o / cto * freq[k]
        bf += math.cos(RD * ttl) ** 2 * freq[k]
        sob += frho * PI / ctscto * freq[k]
        sof += ftau * PI / ctscto * freq[k]
    return ks, ko, bf, sob, sof


def diffuse_scatter_vnir(lai, rho, tau, rs, sza, esun, esky, vza=0.0, raa=0.0, ala=53.0, return_incident=False):
    epsc, epss = 1.0 - rho - tau, 1.0 - rs
    ks, ko, bf, sob, sof = _scattering_coefficients(sza, vza, raa, ala)
    sdb, sdf = 0.5 * (ks + bf), 0.5 * (ks - bf)
    dob, dof = 0.5 * (ko + bf), 0.5 * (ko - bf)
    ddb, ddf = 0.5 * (1.0 + bf), 0.5 * (1.0 - bf)
    sigb, sigf = ddb * rho + ddf * tau, ddf * rho + ddb * tau
    att = 1.0 - sigf
    m = math.sqrt(max((att + sigb) * (att - sigb), 0.0))
    sb, sf = sdb * rho + sdf * tau, sdf * rho + sdb * tau
    vb, vf = dob * rho + dof * tau, dof * rho + dob * tau
    e1, e2 = math.exp(-m * lai), math.exp(-2.0 * m * lai)
    rinf = (att - m) / max(sigb, 1e-12)
    denom = max(1.0 - rinf * rinf * e2, 1e-12)
    re = rinf * e1
    ps = (sf + sb * rinf) * jfunc1(ks, m, lai)
    qs = (sf * rinf + sb) * jfunc2(ks, m, lai)
    pv = (vf + vb * rinf) * jfunc1(ko, m, lai)
    qv = (vf * rinf + vb) * jfunc2(ko, m, lai)
    rho_dd = rinf * (1.0 - e2) / denom
    tau_dd = (1.0 - rinf * rinf) * e1 / denom
    tau_sd = (ps - re * qs) / denom
    rho_sd = (qs - re * ps) / denom
    tau_ss = math.exp(-ks * lai)
    eplu1 = rs * ((tau_ss + tau_sd) * esun + tau_dd * esky) / denom
    e0 = rho_sd * esun + rho_dd * esky + tau_dd * eplu1
    emin1 = tau_sd * esun + tau_dd * esky + rho_dd * eplu1
    leaf_diffuse = max(e0 * epsc, 0.0)
    soil_incident = max(emin1, 0.0)
    soil_diffuse = soil_incident * epss
    if return_incident:
        return leaf_diffuse, soil_diffuse, soil_incident
    return leaf_diffuse, soil_diffuse


def directional_single_scatter(lai, rho, tau, sza, vza, raa, ala=53.0):
    """Directional leaf single-scattering term used by the observation RT pass."""
    lai = max(float(lai), 0.0)
    if lai <= 0.0:
        return 0.0
    ks, ko, _, sob, sof = _scattering_coefficients(float(sza), float(vza), float(raa), ala)
    delta = ks - ko
    if abs(delta) < 1e-8:
        path = lai * math.exp(-ks * lai)
    else:
        path = (math.exp(-ko * lai) - math.exp(-ks * lai)) / delta
    return max((sob * float(rho) + sof * float(tau)) * path, 0.0)


def diffuse_scatter_tir(lai, rho, tau, rs, sza, esky, tss, tsh, tls, tlh, vza=0.0, raa=0.0, ala=53.0):
    epsc, epss = 1.0 - rho - tau, 1.0 - rs
    ks, ko, bf, sob, sof = _scattering_coefficients(sza, vza, raa, ala)
    sdb, sdf = 0.5 * (ks + bf), 0.5 * (ks - bf)
    dob, dof = 0.5 * (ko + bf), 0.5 * (ko - bf)
    ddb, ddf = 0.5 * (1.0 + bf), 0.5 * (1.0 - bf)
    sigb, sigf = ddb * rho + ddf * tau, ddf * rho + ddb * tau
    att = 1.0 - sigf
    m = math.sqrt(max((att + sigb) * (att - sigb), 0.0))
    sb, sf = sdb * rho + sdf * tau, sdf * rho + sdb * tau
    vb, vf = dob * rho + dof * tau, dof * rho + dob * tau
    e1, e2 = math.exp(-m * lai), math.exp(-2.0 * m * lai)
    rinf = (att - m) / max(sigb, 1e-12)
    denom = max(1.0 - rinf * rinf * e2, 1e-12)
    re = rinf * e1
    ps = (sf + sb * rinf) * jfunc1(ks, m, lai)
    qs = (sf * rinf + sb) * jfunc2(ks, m, lai)
    pv = (vf + vb * rinf) * jfunc1(ko, m, lai)
    qv = (vf * rinf + vb) * jfunc2(ko, m, lai)
    tau_dd = (1.0 - rinf * rinf) * e1 / denom
    tau_sd = (ps - re * qs) / denom
    rho_dd = rinf * (1.0 - e2) / denom
    f_hs = (1.0 - rinf * rinf) * (1.0 - rs) / max(1.0 - rinf * rs, 1e-12)
    f_hc = m * (1.0 - rinf)
    f_bottom = (rs - rinf) / max(1.0 - rinf * rs, 1e-12)
    es = 0.5 * (SIGMA * tss**4 + SIGMA * tsh**4)
    ev = 0.5 * (SIGMA * tls**4 + SIGMA * tlh**4)
    f11 = esky * (1.0 - f_hc) + ev * f_hc
    f21 = f_bottom * f11 + es * f_hs
    emin = (f11 + rinf * f21) / max(1.0 - rinf * rinf, 1e-12)
    eplu = (f21 + rinf * f11) / max(1.0 - rinf * rinf, 1e-12)
    return max(emin * epss, 0.0), max((eplu + esky) * epsc, 0.0)


def tirt_direct(lai, hspot, sza, saa, vza, vaa):
    psoil = float(np.asarray(gap_probability_hom_analytical(lai, vza)).reshape(-1)[0])
    pleaf = 1.0 - psoil
    vsa = np.asarray([abs(vaa - saa)])
    soil_sun = float(np.asarray(hotspot_layer(lai, hspot, np.asarray([vza]), np.asarray([sza]), vsa)).reshape(-1)[0])
    leaf_sun = float(np.asarray(hotspot_analytical(lai, hspot, np.asarray([vza]), np.asarray([sza]), vsa)).reshape(-1)[0])
    return soil_sun, psoil - soil_sun, leaf_sun, pleaf - leaf_sun


def tirt_scatter(lai, rho_soil, rho_leaf, sza):
    bv = float(np.asarray(gap_probability_hom_analytical(lai, 0.0)).reshape(-1)[0])
    m = float(np.asarray(gap_probability_hom_analytical(lai, 0.0)).reshape(-1)[0])
    m = math.exp(-0.825 * lai * 0.5)
    alpha = [0.2855375, 0.2885375, 0.2964427, 0.3003953, 0.3083004, 0.3201581, 0.3399209, 0.3715415, 0.4189723, 1.0]
    ms = (1.0 - rho_leaf) * (bv * (1.0 - m) * rho_soil + (1.0 - alpha[0]) * (1.0 - bv * m) * (1.0 - bv) * rho_leaf)
    volume = float(np.asarray(hotspot_vegetation_volume(lai, np.asarray([sza]))).reshape(-1)[0])
    return 0.0, 0.0, max(ms * volume, 0.0), max(ms * (1.0 - volume), 0.0)
