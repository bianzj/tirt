"""Closed-form four-component fractions of Poisson-distributed turbid spheroid crowns.

Soil view: exact Poisson-spheroid blocking, J = A [1 - 2(1 - (1+x)e^-x)/x^2], x = 1.5 k l_mean.
Sunlit soil: the sun blocking not already shared with the view ray; the shared part is the lens
overlap of the two crown projections, whose centres are H_c |tan(v) - tan(s)| apart (H_c = crown
base + vertical semiaxis), times the probability that one crown blocks both rays.
Sunlit leaves: lit share of the visible crown surface from the phase angle in sphere-normalised
space, sunlight through the crown for the far side, and the finite depth of the visible layer.
All angles are evaluated together as arrays.
"""
import numpy as np

G = 0.5
# Fitted once to the exact Poisson-crown integral; not tuned to HiSTREAM.
LEAF_COEFFICIENTS = (0.567, 0.369, 1.595)


def _blocking(k, rcr, hcr, theta):
    area = np.pi * rcr * np.sqrt(rcr ** 2 + (hcr * np.tan(theta)) ** 2)
    volume = 4.0 / 3.0 * np.pi * rcr ** 2 * hcr
    depth = k * volume / (area * np.cos(theta))
    x = np.maximum(1.5 * depth, 1.0e-12)
    safe = np.where(x > 1.0e-4, x, 1.0)
    mean_t = np.where(x > 1.0e-4, 2.0 * (1.0 - (1.0 + safe) * np.exp(-safe)) / safe ** 2, 1.0 - x * 2.0 / 3.0)
    return area * (1.0 - mean_t), area, mean_t, depth


def _lens(r1, r2, d):
    small, big = np.minimum(r1, r2), np.maximum(r1, r2)
    inside = d <= big - small
    apart = d >= big + small
    dd = np.where(inside | apart, 1.0, d)
    a1 = small ** 2 * np.arccos(np.clip((dd * dd + small * small - big * big) / (2 * dd * small), -1, 1))
    a2 = big ** 2 * np.arccos(np.clip((dd * dd + big * big - small * small) / (2 * dd * big), -1, 1))
    a3 = 0.5 * np.sqrt(np.maximum((-dd + small + big) * (dd + small - big) * (dd - small + big) * (dd + small + big), 0.0))
    return np.where(inside, np.pi * small ** 2, np.where(apart, 0.0, a1 + a2 - a3))


def crown_fractions(lai, std, hspot, hcr, rcr, hcb, vza, sza, raa):
    """(4, N) sunlit soil, shaded soil, sunlit leaf, shaded leaf for arrays of angles (degrees)."""
    vza, sza, raa = (np.asarray(a, float).reshape(-1) for a in np.broadcast_arrays(vza, sza, raa))
    if lai <= 0.0 or std <= 0.0:
        return np.vstack([np.ones_like(vza), np.zeros_like(vza), np.zeros_like(vza), np.zeros_like(vza)])
    rd = np.pi / 180.0
    tv, ts = np.minimum(np.abs(vza), 89.9) * rd, np.minimum(np.abs(sza), 89.9) * rd
    phi = np.abs((raa + 180.0) % 360.0 - 180.0) * rd
    k = G * lai / (std * 4.0 / 3.0 * np.pi * rcr ** 2 * hcr)
    j_v, a_v, t_v, tau_v = _blocking(k, rcr, hcr, tv)
    j_s, a_s, t_s, tau_s = _blocking(k, rcr, hcr, ts)

    centre = max(hcb, 0.0) + hcr
    tan_v, tan_s = np.tan(tv), np.tan(ts)
    dtan = np.sqrt(np.maximum(tan_v ** 2 + tan_s ** 2 - 2 * tan_v * tan_s * np.cos(phi), 0.0))
    cv, cs = np.cos(tv), np.cos(ts)
    a = dtan * 2.0 / (G / cv + G / cs) / max(hspot * 2.0 * hcr, 1.0e-9) * 4.0 / 3.0 * hcr
    w = np.where(a > 1.0e-9, (1.0 - np.exp(-a)) / np.where(a > 1.0e-9, a, 1.0), 1.0)
    # joint crown transmittance, evaluated in log space; it can never exceed either single transmittance
    log_t_vs = (np.log(np.maximum(t_v, 1.0e-300)) + np.log(np.maximum(t_s, 1.0e-300))
                + np.minimum(np.minimum(np.sqrt(tau_v * tau_s) * w, tau_v), tau_s))
    t_vs = np.minimum(np.minimum(np.exp(np.minimum(log_t_vs, 0.0)), t_v), t_s)
    shared = _lens(np.sqrt(a_v / np.pi), np.sqrt(a_s / np.pi), centre * dtan) * np.maximum(1.0 - t_v - t_s + t_vs, 0.0)
    shared = np.minimum(np.minimum(shared, j_v), j_s)
    soil = np.exp(-std * j_v)
    soil_sunlit = soil * np.exp(-std * (j_s - shared))

    c_depth, c_dark, c_lit = LEAF_COEFFICIENTS
    # phase angle after mapping the spheroid onto a sphere (z scaled by rcr/hcr)
    q = rcr / hcr
    v = np.stack([np.sin(tv) * np.cos(phi), np.sin(tv) * np.sin(phi), np.cos(tv) * q])
    s = np.stack([np.sin(ts), np.zeros_like(ts), np.cos(ts) * q])
    cos_phase = np.sum(v * s, 0) / np.linalg.norm(v, axis=0) / np.linalg.norm(s, axis=0)
    lit = (0.5 * (1.0 + cos_phase)) ** c_lit
    r_v = k / np.sqrt(np.sin(tv) ** 2 / rcr ** 2 + np.cos(tv) ** 2 / hcr ** 2)
    r_s = k / np.sqrt(np.sin(ts) ** 2 / rcr ** 2 + np.cos(ts) ** 2 / hcr ** 2)
    ratio = c_depth * (1.0 - w) * r_s / np.maximum(r_v, 1.0e-9)
    tau = np.maximum(2.0 * r_v, 1.0e-9)
    depth = (1.0 - np.exp(-tau * (1.0 + ratio))) / (1.0 + ratio) / (1.0 - np.exp(-tau))
    share = np.minimum((lit + (1.0 - lit) * np.exp(-c_dark * r_s * (1.0 - w))) * depth, 1.0)
    leaf = 1.0 - soil
    return np.vstack([soil_sunlit, soil - soil_sunlit, leaf * share, leaf * (1.0 - share)])
