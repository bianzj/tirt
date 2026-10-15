"""Periodic row optics: condition the first collision on open corridors.

The diffuse closure stays separate from the joint sun/view geometry. No
thermal row probabilities or thermal scattering functions are modified.
"""
import math

from .vnir import spherical_coefficients, joint_transmission


def clear_corridor(width, blank, height, sza, saa, vza, vaa, row_azimuth):
    """Fraction traversing neither foliage nor its projected sun/view shadow."""
    if width <= 0 or blank < 0 or height <= 0:
        raise ValueError('Invalid row dimensions.')
    sun = height*math.tan(math.radians(sza))*math.sin(math.radians(saa-row_azimuth))
    view = height*math.tan(math.radians(vza))*math.sin(math.radians(vaa-row_azimuth))
    shadow = max(0., sun, view)-min(0., sun, view)
    return max(blank-shadow, 0.)/(width+blank)


def conditioned_first(lai, hspot, rho, tau, soil, sza, vza, raa,
                      solar_gap, view_gap, clear):
    """Joint first scattering with the empty phase removed before extinction.

    The remaining projected row domain is an effective slab. Along a row,
    this reduces exactly to clear soil plus a slab with LAI/row-cover. The
    caller retains the established closure where no empty corridor survives,
    and blends the two closures with the projected corridor survival fraction.
    """
    ks, ko, sob, sof = spherical_coefficients(float(sza), float(vza), float(raa))
    cover = 1-clear
    gs = min(1., max((solar_gap-clear)/cover, 1e-12))
    gv = min(1., max((view_gap-clear)/cover, 1e-12))
    depth = max(-math.log(gs)/ks, 1e-10)
    extinction = -math.log(gv)/depth
    joint, path = joint_transmission(depth, ks, extinction, hspot, sza, vza, raa)
    leaf = cover*(sob*rho+sof*tau)*extinction/ko*path
    ground = soil*(clear+cover*joint)
    return leaf+ground


def primary_correction(settings, sza, saa, vza, vaa, rho, tau, soil, solar_gap, fractions):
    clear = clear_corridor(*[float(settings[key]) for key in ('row_width','row_blank','row_height')],
                           sza, saa, vza, vaa, float(settings['row_azimuth']))
    if clear <= 1e-14 or float(settings['lai']) <= 1e-10:
        return 0.
    view_gap = float(fractions[0]+fractions[1])
    raa = abs((vaa-saa+180.) % 360.-180.)
    args = (float(settings['lai']), float(settings['hspot']), rho, tau, soil,
            sza, vza, raa, solar_gap, view_gap)
    _, ko, sob, sof = spherical_coefficients(float(sza), float(vza), float(raa))
    old = (sob*rho+sof*tau)/ko*fractions[2]+soil*fractions[0]
    # Projected corridor survival is the interpolation weight. It is one
    # when both rays follow the rows and zero when no corridor survives.
    # This keeps the established closed-row approximation continuous while
    # recovering the conditional slab exactly in the open, along-row limit.
    weight = clear*(float(settings['row_width'])+float(settings['row_blank']))/float(settings['row_blank'])
    return weight*(conditioned_first(*args, clear)-old)
