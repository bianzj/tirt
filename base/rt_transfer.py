import math


PI = math.pi
RD = PI / 180.0


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
