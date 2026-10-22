"""Closed-form four-component fractions of row crops (turbid rows of width w, gap b, height H).

In the cross-row plane a ray crossing the canopy layer moves D = H tan(theta) sin(phi - row) sideways;
its in-row path is (H/|D|) c(x), where c(x) is the row length covered by the window [x, x+D]. c(x) is
piecewise linear in the start position x, so the row-period means of exp(-tau) are exact sums over the
breakpoints. Sunlit visible leaves use the same exact x-means for leaves at each depth below the row
top, with a fixed 8-point Gauss rule over depth.
"""
import numpy as np

G = 0.5


def _covered(x, length, period, width):
    """Row length inside [x, x + length] (length may be negative)."""
    def cumulative(t):
        return np.floor(t / period) * width + np.minimum(np.mod(t, period), width)
    return np.abs(cumulative(x + length) - cumulative(x))


def _breakpoints(shift, period, width):
    """Start positions in [0, period) where a window edge crosses a row edge."""
    edges = np.array([0.0, width])
    return np.mod(np.concatenate([edges, edges - shift]), period)


def _optical_depth(x, shift, height, extinction, mu, period, width):
    if abs(shift) < 1.0e-9:
        return extinction / mu * height * (np.mod(x, period) < width)
    return extinction / mu * height / abs(shift) * _covered(x, shift, period, width)


def _period_mean_exp(exponent, points, period, upper=None):
    """Exact mean over one period of exp(e(x)) for e linear between consecutive points.

    With `upper`, only x in [0, upper) contributes (still divided by the period).
    """
    end = period if upper is None else upper
    x = np.unique(np.concatenate([[0.0, end], points[(points > 0.0) & (points < end)]]))
    mid = 0.5 * (x[:-1] + x[1:])
    # evaluate each piece from the inside to avoid the jump of a discontinuous piece at its ends
    eps = 1.0e-9 * period
    left, right = exponent(x[:-1] + eps), exponent(x[1:] - eps)
    span = x[1:] - x[:-1]
    diff = right - left
    piece = np.where(np.abs(diff) > 1.0e-10, (np.exp(right) - np.exp(left)) / np.where(np.abs(diff) > 1.0e-10, diff, 1.0),
                     np.exp(exponent(mid)))
    return float(np.sum(span * piece) / period)


_HEIGHT_NODES, _HEIGHT_WEIGHTS = np.polynomial.legendre.leggauss(8)
_HEIGHT_NODES, _HEIGHT_WEIGHTS = 0.5 * (_HEIGHT_NODES + 1.0), 0.5 * _HEIGHT_WEIGHTS


def row_fractions(lai, hspot, row_width, row_blank, row_height, vza, vaa, sza, saa, row_azimuth):
    """Sunlit soil, shaded soil, sunlit leaf, shaded leaf for one direction (degrees)."""
    if lai <= 0.0:
        return np.array([1.0, 0.0, 0.0, 0.0])
    rd = np.pi / 180.0
    tv, ts = min(abs(vza), 89.9) * rd, min(abs(sza), 89.9) * rd
    cv, cs = np.cos(tv), np.cos(ts)
    period = row_width + row_blank
    extinction = G * lai * period / (row_width * row_height)
    shift_v = row_height * np.tan(tv) * np.sin((vaa - row_azimuth) * rd)
    shift_s = row_height * np.tan(ts) * np.sin((saa - row_azimuth) * rd)
    points = np.concatenate([_breakpoints(shift_v, period, row_width), _breakpoints(shift_s, period, row_width)])

    def tau_v(x):
        return _optical_depth(x, shift_v, row_height, extinction, cv, period, row_width)

    def tau_s(x):
        return _optical_depth(x, shift_s, row_height, extinction, cs, period, row_width)

    # leaf-scale hotspot overlap, same weight as hotspot_layer
    delta = np.tan(ts) ** 2 + np.tan(tv) ** 2 - 2 * np.tan(ts) * np.tan(tv) * np.cos((saa - vaa) * rd)
    alpha = np.sqrt(max(delta, 0.0)) * 2.0 / (G / cv + G / cs) / max(hspot, 1.0e-12)
    w = (1.0 - np.exp(-alpha)) / alpha if alpha > 1.0e-9 else 1.0

    soil = _period_mean_exp(lambda x: -tau_v(x), points, period)
    soil_sunlit = _period_mean_exp(lambda x: -tau_v(x) - tau_s(x) + np.sqrt(tau_v(x) * tau_s(x)) * w, points, period)
    soil_sunlit = min(soil_sunlit, soil)

    # Leaves at depth h below the row top: x-mean over the row strip, 8-point Gauss in h.
    leaf_sunlit = 0.0
    for node, weight in zip(_HEIGHT_NODES, _HEIGHT_WEIGHTS):
        depth = node * row_height
        sv, ss = shift_v * node, shift_s * node
        level = np.concatenate([_breakpoints(sv, period, row_width), _breakpoints(ss, period, row_width)])

        def exponent(x, sv=sv, ss=ss, depth=depth):
            a = _optical_depth(x, sv, depth, extinction, cv, period, row_width)
            b = _optical_depth(x, ss, depth, extinction, cs, period, row_width)
            return -a - b + np.sqrt(a * b) * w
        leaf_sunlit += weight * row_height * extinction / cv * _period_mean_exp(exponent, level, period, upper=row_width)
    leaf = 1.0 - soil
    leaf_sunlit = min(leaf_sunlit, leaf)
    return np.array([soil_sunlit, soil - soil_sunlit, leaf_sunlit, leaf - leaf_sunlit])
