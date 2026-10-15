"""Diffuse eigenmode source profile for whole-crown optical component exchange."""
import math
from functools import lru_cache

from scipy.optimize import brentq

from .gap import gap_probability_crown_analytical
from .optical_hotspot import attenuation_mean


def escape_ratios(view, depth, attenuation):
    view_extinction = .5/math.cos(math.radians(view))
    normalization = attenuation_mean(attenuation*depth)*attenuation_mean(view_extinction*depth)
    upward = attenuation_mean((attenuation+view_extinction)*depth)/normalization
    numerator = (math.exp(-min(attenuation, view_extinction)*depth)*
                 attenuation_mean(abs(view_extinction-attenuation)*depth))
    return upward, numerator/normalization


@lru_cache(maxsize=512)
def beam_down_escape(shape, solar, order):
    """Downward escape of a beam-born source, rather than a uniform source.

    Use the same equivalent crown depth and normalized directional escape
    as the higher-order profile, with direct extinction G/cos(solar). No
    material reflectance or reference BRF enters this first-source geometry.
    """
    from .optical_crown import quadrature
    depth = source_profile(shape, 0., 0., order)[0]
    attenuation = .5/math.cos(math.radians(solar))
    cosine, weights, _ = quadrature(order)
    escape = 0.
    for value, weight in zip(cosine, weights):
        view = math.degrees(math.acos(value))
        gap = float(gap_probability_crown_analytical(*shape[:4], view))
        escape += 2*value*weight*(1-gap)/shape[0]*escape_ratios(view, depth, attenuation)[1]
    if not 0 <= escape <= 1:
        raise ValueError('Invalid beam-born downward escape probability.')
    return float(escape)


@lru_cache(maxsize=512)
def source_profile(shape, rho, tau, order):
    from .optical_crown import geometry, quadrature
    escape, _ = geometry(shape, order)
    depth = (1e-10 if escape >= 1-1e-12 else
             brentq(lambda value: attenuation_mean(value)-escape, 1e-10, max(4., 2/escape)))
    backward, forward = (2*rho+tau)/3, (rho+2*tau)/3
    attenuation = math.sqrt(max(0., (1-forward)**2-backward**2))
    cosine, weights, _ = quadrature(order)
    upward_escape, downward_escape = 0., 0.
    for value, weight in zip(cosine, weights):
        zenith = math.degrees(math.acos(value))
        gap = float(gap_probability_crown_analytical(*shape[:4], zenith))
        upward, downward = escape_ratios(zenith, depth, attenuation)
        probability = 2*value*weight*(1-gap)/shape[0]
        upward_escape += probability*upward
        downward_escape += probability*downward
    if not (0 <= upward_escape <= 1+1e-10 and 0 <= downward_escape <= 1+1e-10):
        raise ValueError('Diffuse crown source produces invalid escape probabilities.')
    return depth, attenuation, upward_escape, downward_escape
