"""Node radiosity for multiple scattering with closed-form, reciprocal exchange factors.

Nodes carry a capacity A (soil: ground fraction; leaves: 4 G LAI per unit ground), so that
A_i F_ij = A_j F_ji. With X(tau) = 1 - 2 E3(tau) the diffuse interception of a turbid layer of
optical depth tau = G LAI, every factor follows from reciprocity:
  layer -> sky (top) = X(tau)/A, layer self = 1 - X(tau)/(2 tau),
  adjacent layers U over D: A_U F_UD = X(tau_U) + X(tau_D) - X(tau_U + tau_D).
Schedules:
  hom2   : soil + one canopy layer (p-theory with soil coupling, p = 1 - X(tau)/(2 tau))
  hom3   : soil + lower + upper canopy layers
  row2/3 : soil (or covered + open soil) + canopy for rows, open ground = gap / period
  crown2/3 : the same for Poisson crowns, open ground = exp(-density pi a^2)
For the split schedules the CI two-phase surrogate puts every clear path on open ground, so
the covered and open soil gaps average back to the CI hemispherical gap exactly.
Sources, interception and view weights per node use the same CI/Beer-Lambert gaps.
"""
import functools
import math

import numpy as np
from scipy.special import expn

from . import row_closure as rc

G = 0.5
_HEMI_MU, _HEMI_WEIGHT = np.polynomial.legendre.leggauss(32)
_HEMI_MU, _HEMI_WEIGHT = (_HEMI_MU + 1.) / 2., _HEMI_WEIGHT / 2.
# Row reflection symmetry reduces full 96-azimuth integration to one quadrant.
_HEMI_AZIMUTHS = (np.arange(24) + .5) * (90. / 24.)


def _intercept(tau):
    return 1.0 - 2.0 * float(expn(3, tau)) if tau > 0.0 else 0.0


# Inverting a gap of 1 - O(LAI) for the CI clear fraction loses all precision as LAI -> 0,
# while that fraction has a finite LAI -> 0 limit; keep its value at this LAI below it.
_PHASE_LAI = 1.0e-4


def _phases(lai, gap_of, zenith, relative, mu):
    """CI clear fraction and covered optical depth along one direction."""
    if lai >= _PHASE_LAI:
        return rc._phases(gap_of(lai)(zenith, relative), G * lai / mu)
    clear, _ = rc._phases(gap_of(_PHASE_LAI)(zenith, relative), G * _PHASE_LAI / mu)
    return (clear, G * lai / mu / (1.0 - clear)) if clear < 1.0 else (1.0, 0.0)


class Nodes:
    """Capacities, exchange factors, soil/leaf type and height level of each node.

    Scattering is directional: soil reflects upward; a canopy node keeps an up- and a down-going
    state. Leaves send back = (2 rho + tau)/3 to the incident side and forward = (rho + 2 tau)/3 through
    (spherical bi-Lambertian, as in the four-stream closure). An emitting canopy node sends half of its
    power each way; exchange toward higher levels and the sky comes from the up half, toward lower
    levels and the soil from the down half, and each half's remainder is self-interception.
    """

    def __init__(self, capacity, factors, soil, level):
        self.capacity, self.factors = np.asarray(capacity, float), np.asarray(factors, float)
        self.soil, self.level = np.asarray(soil, bool), np.asarray(level, float)
        if np.any(self.factors.sum(axis=1) > 1.0 + 1e-9):
            raise ValueError("Node exchange factors exceed unity.")
        n = len(self.capacity)
        # States: every node's up state, then the down state of each canopy node.
        self.up = np.arange(n)
        self.down = n + np.cumsum(~self.soil) - 1
        self.size = n + int((~self.soil).sum())
        # Arrivals: node * 2 + side, side 0 = from below, 1 = from above.
        self.arrive = np.zeros((self.size, 2 * n))
        escape = self.escape()
        for i in range(n):
            higher, lower = self.level > self.level[i], self.level < self.level[i]
            if self.soil[i]:
                self.arrive[self.up[i], 2 * np.arange(n)] = self.factors[i]
                continue
            up, down = self.up[i], self.down[i]
            self.arrive[up, 2 * np.flatnonzero(higher)] = 2.0 * self.factors[i, higher]
            self.arrive[up, 2 * i] = max(0.0, 1.0 - 2.0 * (self.factors[i, higher].sum() + escape[i]))
            self.arrive[down, 2 * np.flatnonzero(lower) + 1] = 2.0 * self.factors[i, lower]
            self.arrive[down, 2 * i + 1] = max(0.0, 1.0 - 2.0 * self.factors[i, lower].sum())
        self.radiance = np.zeros(self.size)
        self.radiance[self.up] = np.where(self.soil, 1.0, 2.0) / (np.pi * self.capacity)

    def scattering(self, soil_albedo, rho, tau):
        """Arrivals -> leaving power of each state."""
        return self._scatter(soil_albedo, (2.0 * rho + tau) / 3.0, (rho + 2.0 * tau) / 3.0)

    def _scatter(self, soil_albedo, back, forward):
        n = len(self.capacity)
        out = np.zeros((self.size, 2 * n))
        for i in range(n):
            if self.soil[i]:
                out[self.up[i], 2 * i:2 * i + 2] = soil_albedo
            else:
                out[self.up[i], 2 * i], out[self.up[i], 2 * i + 1] = forward, back
                out[self.down[i], 2 * i], out[self.down[i], 2 * i + 1] = back, forward
        return out

    def solve(self, soil_albedo, rho, tau, first):
        """All orders of leaving power from first-leaving power `first` (states x k)."""
        transfer = self.scattering(soil_albedo, rho, tau) @ self.arrive.T
        return np.linalg.solve(np.eye(self.size) - transfer, first)

    def emission(self, per_capacity):
        """Leaving power of an emission given per unit capacity (nodes x k); canopy splits up/down."""
        power = np.pi * self.capacity[:, None] * per_capacity
        states = np.zeros((self.size, power.shape[1]))
        states[self.up] = np.where(self.soil[:, None], power, 0.5 * power)
        states[self.down[~self.soil]] = 0.5 * power[~self.soil]
        return states

    def escape(self):
        return np.clip(1.0 - self.factors.sum(axis=1), 0.0, 1.0)


@functools.lru_cache(maxsize=64)
def homogeneous(lai, layers):
    """Soil + `layers` equal canopy layers (top layer last)."""
    tau = G * lai / layers
    capacity = [1.0] + [4.0 * tau] * layers
    n = layers + 1
    exchange = np.zeros((n, n))
    for i in range(1, n):
        for j in range(1, n):
            if i == j:
                exchange[i, j] = 4.0 * tau - 2.0 * _intercept(tau)
            elif abs(i - j) == 1:
                exchange[i, j] = 2.0 * _intercept(tau) - _intercept(2.0 * tau)
            else:
                gap = abs(i - j) - 1
                exchange[i, j] = (_intercept((gap + 1) * tau) - _intercept(gap * tau)) - \
                                 (_intercept((gap + 2) * tau) - _intercept((gap + 1) * tau))
        below = (i - 1) * tau
        exchange[0, i] = exchange[i, 0] = _intercept(below + tau) - _intercept(below)
    return Nodes(capacity, exchange / np.asarray(capacity)[:, None], [True] + [False] * layers, np.arange(n))


def _soil_gaps(lai, gap_of, open_fraction):
    """Cosine-weighted hemispherical gaps of covered and open soil from the CI two-phase surrogate.

    Every clear path starts on open ground, so open_fraction * M_open + (1 - open_fraction) * M_covered = M.
    """
    covered = opened = 0.0
    zenith = np.degrees(np.arccos(_HEMI_MU))
    for mu, weight, z in zip(_HEMI_MU, _HEMI_WEIGHT, zenith):
        u = g = 0.0
        for azimuth in _HEMI_AZIMUTHS:
            clear, depth = _phases(lai, gap_of, z, azimuth, mu)
            share = min(1.0, clear / open_fraction)
            u += math.exp(-depth)
            g += share + (1.0 - share) * math.exp(-depth)
        covered += 2.0 * weight * mu * u / len(_HEMI_AZIMUTHS)
        opened += 2.0 * weight * mu * g / len(_HEMI_AZIMUTHS)
    return covered, opened


def _covered_nodes(lai, gap_of, open_fraction, split_soil):
    """Soil (covered + open, or one node) + canopy, with reciprocal closed-form exchange factors."""
    canopy = 4.0 * G * lai
    if split_soil and 0.0 < open_fraction < 1.0:
        covered, opened = _soil_gaps(lai, gap_of, open_fraction)
        area = np.array([1.0 - open_fraction, open_fraction])
        to_canopy = np.array([1.0 - covered, 1.0 - opened])
    else:
        gap, hemi = gap_of(lai), 0.0
        for mu, weight, z in zip(_HEMI_MU, _HEMI_WEIGHT, np.degrees(np.arccos(_HEMI_MU))):
            hemi += 2.0 * weight * mu * np.mean([gap(z, a) for a in _HEMI_AZIMUTHS])
        area, to_canopy = np.array([1.0]), np.array([1.0 - hemi])
    intercepted = float(np.sum(area * to_canopy))
    n = len(area) + 1
    factors = np.zeros((n, n))
    factors[:-1, -1] = to_canopy
    factors[-1, :-1] = area * to_canopy / canopy
    factors[-1, -1] = 1.0 - intercepted / (2.0 * G * lai)
    return Nodes(np.r_[area, canopy], factors, [True] * len(area) + [False], [0] * len(area) + [1])


def _covered_split(lai, gap_of, open_fraction, split_soil, zenith, relative):
    """Per-node fractions along one direction: soil nodes (gap to ground) and the canopy (collided)."""
    p = gap_of(lai)(zenith, relative)
    if not (split_soil and 0.0 < open_fraction < 1.0):
        return np.array([p, 1.0 - p])
    clear, depth = _phases(lai, gap_of, zenith, relative, math.cos(math.radians(zenith)))
    share = min(1.0, clear / open_fraction)
    covered = (1.0 - open_fraction) * math.exp(-depth)
    opened = open_fraction * (share + (1.0 - share) * math.exp(-depth))
    # At grazing incidence both transmissions can underflow together. All of the
    # beam then strikes canopy; avoid the otherwise indeterminate 0/0 soil split.
    soil = np.array([covered, opened])
    total = covered + opened
    soil = soil * p / total if total > 0.0 else np.zeros(2)
    return np.r_[soil, 1.0 - p]


def _row_gap(lai, width, blank, height):
    return lambda z, a: rc._gap(lai, width, blank, height, z, a)


@functools.lru_cache(maxsize=64)
def row(lai, width, blank, height, split_soil):
    return _covered_nodes(lai, lambda x: _row_gap(x, width, blank, height), blank / (width + blank), split_soil)


def _crown_gap(lai, density, vertical, horizontal):
    """Directional gap of the Poisson-spheroid crowns, the same gap as the crown four components."""
    from .crown_parametric import crown_fractions

    def gap(z, a):
        parts = crown_fractions(lai, density, 0.1, vertical, horizontal, vertical, z, z, 0.0)[:, 0]
        return float(parts[0] + parts[1])
    return gap


def crown_open_fraction(density, horizontal):
    """Ground outside the vertical projections of Poisson crowns."""
    return math.exp(-density * math.pi * horizontal * horizontal)


@functools.lru_cache(maxsize=64)
def crown(lai, density, vertical, horizontal, split_soil):
    return _covered_nodes(lai, lambda x: _crown_gap(x, density, vertical, horizontal),
                          crown_open_fraction(density, horizontal), split_soil)


def _hom_split(lai, layers, zenith):
    """Per-node fractions along one direction for soil + equal layers (top last)."""
    k = G / math.cos(math.radians(zenith))
    tops = np.exp(-k * lai * np.arange(layers - 1, -1, -1) / layers)  # transmission down to each layer top
    hits = tops * (1.0 - math.exp(-k * lai / layers))
    return np.r_[math.exp(-k * lai), hits]


def _scheme(kind, lai, shape):
    """shape: rows (width, blank, height); crowns (density, vertical semiaxis, horizontal radius)."""
    if kind == "hom2":
        return homogeneous(lai, 1), lambda z, r: _hom_split(lai, 1, z)
    if kind == "hom3":
        return homogeneous(lai, 2), lambda z, r: _hom_split(lai, 2, z)
    split = kind in ("row3", "crown3")
    if kind in ("crown2", "crown3"):
        density, vertical, horizontal = shape
        gap_of, open_fraction = (lambda x: _crown_gap(x, *shape)), crown_open_fraction(density, horizontal)
        return crown(lai, *shape, split), lambda z, r: _covered_split(lai, gap_of, open_fraction, split, z, r)
    width, blank, height = shape
    gap_of, open_fraction = (lambda x: _row_gap(x, *shape)), blank / (width + blank)
    return row(lai, *shape, split), lambda z, r: _covered_split(lai, gap_of, open_fraction, split, z, r)


def thermal_scatter(kind, lai, shape, emis_leaf, emis_soil, vza, vaa, sza, saa, row_azimuth=0.0):
    """Multiple-scattering effective emissivities [soil sun, soil shade, leaf sun, leaf shade]."""
    if lai <= 0.0:
        return np.zeros(4)
    nodes, split = _scheme(kind, float(lai), tuple(shape))
    zv, zs = rc._zenith(vza), rc._zenith(sza)
    view = split(zv, float(vaa) - row_azimuth)
    sun = split(zs, float(saa) - row_azimuth)
    area = np.where(nodes.soil, nodes.capacity, 0.0)
    leaf_area = np.where(nodes.soil, 0.0, nodes.capacity / (4.0 * G))
    # Sunlit share: soil = beam reaching the node per node area; leaves = intercepted beam per leaf area.
    k_sun = G / math.cos(math.radians(zs))
    sunlit = np.where(nodes.soil, np.divide(sun, area, out=np.zeros_like(sun), where=area > 0.0),
                      np.divide(sun, k_sun * leaf_area, out=np.zeros_like(sun), where=leaf_area > 0.0))
    sunlit = np.clip(sunlit, 0.0, 1.0)
    source = np.zeros((len(sunlit), 4))
    source[nodes.soil, 0], source[nodes.soil, 1] = emis_soil * sunlit[nodes.soil], emis_soil * (1.0 - sunlit[nodes.soil])
    source[~nodes.soil, 2], source[~nodes.soil, 3] = emis_leaf * sunlit[~nodes.soil], emis_leaf * (1.0 - sunlit[~nodes.soil])
    emitted = nodes.emission(source)
    scattered = nodes.solve(1.0 - emis_soil, 1.0 - emis_leaf, 0.0, emitted) - emitted
    return _seen(nodes, view) @ scattered


def _seen(nodes, view):
    """Sensor weight of each state: view fraction of the node times its upward radiance per power."""
    weights = np.zeros(nodes.size)
    weights[nodes.up] = view
    return weights * nodes.radiance


def _from_above(nodes, power):
    arrivals = np.zeros(2 * len(nodes.capacity))
    arrivals[2 * np.arange(len(nodes.capacity)) + 1] = power
    return arrivals


def _escaping(nodes):
    """Fraction of each state's leaving power that escapes to the sky."""
    weights = np.zeros(nodes.size)
    weights[nodes.up] = np.where(nodes.soil, 1.0, 2.0) * nodes.escape()
    return weights


def _seen_split(nodes, seen):
    """Sensor weights split into [canopy nodes, soil nodes]."""
    canopy, soil = seen.copy(), seen.copy()
    canopy[nodes.up[nodes.soil]] = 0.0
    soil[nodes.up[~nodes.soil]] = 0.0
    soil[nodes.down[~nodes.soil]] = 0.0
    return np.vstack([canopy, soil])


def beam_response(kind, lai, shape, rho, tau, soil_rho, vza, vaa, sza, saa, row_azimuth=0.0,
                  canopy_weight=1.0, soil_weight=1.0):
    """Unit beam: orders >= 2 toward the view as BRF [through canopy, through soil], and escaped albedo (all orders).

    canopy_weight / soil_weight scale the beam reaching the canopy and the ground (e.g. building shadow).
    """
    nodes, split = _scheme(kind, float(lai), tuple(shape))
    zv, zs = rc._zenith(vza), rc._zenith(sza)
    seen = _seen_split(nodes, _seen(nodes, split(zv, float(vaa) - row_azimuth)))
    weights = np.where(nodes.soil, soil_weight, canopy_weight)
    first = nodes.scattering(soil_rho, rho, tau) @ _from_above(nodes, weights * split(zs, float(saa) - row_azimuth))
    total = nodes.solve(soil_rho, rho, tau, first[:, None])[:, 0]
    return np.pi * (seen @ (total - first)), float(_escaping(nodes) @ total)


def sky_response(kind, lai, shape, rho, tau, soil_rho, vza, vaa, row_azimuth=0.0):
    """Unit isotropic sky, all orders: HDRF toward the view [through canopy, through soil] and diffuse albedo."""
    nodes, split = _scheme(kind, float(lai), tuple(shape))
    seen = _seen_split(nodes, _seen(nodes, split(rc._zenith(vza), float(vaa) - row_azimuth)))
    sky = nodes.escape() * nodes.capacity
    first = nodes.scattering(soil_rho, rho, tau) @ _from_above(nodes, sky / sky.sum())
    total = nodes.solve(soil_rho, rho, tau, first[:, None])[:, 0]
    return np.pi * (seen @ total), float(_escaping(nodes) @ total)


def canopy_nodes(kind, lai, shape, sza, saa, row_azimuth=0.0):
    """Nodes of a schedule and the unit horizontal beam reaching each node first (per unit ground)."""
    nodes, split = _scheme(kind, float(lai), tuple(shape))
    return nodes, split(rc._zenith(sza), float(saa) - row_azimuth)


def from_above(nodes, power):
    """Arrivals of a downward source distributed over the nodes as `power` (per unit ground)."""
    return _from_above(nodes, np.asarray(power, float))


def sky_arrivals(nodes):
    """Arrivals of a unit isotropic sky irradiance (per unit ground), by reciprocity with the sky view."""
    sky = nodes.escape() * nodes.capacity
    return _from_above(nodes, sky / sky.sum())


def absorbed(nodes, soil_albedo, rho, tau, arrivals, sources=None):
    """Power absorbed by each node (bands x nodes): (from the external arrivals directly, from all scattering).

    soil_albedo, rho and tau are scalars or per-band arrays; arrivals (2 n) or (bands, 2 n); sources are leaving
    powers (states) or (bands, states), e.g. Nodes.emission, which also include self-interception.
    """
    soil_albedo, rho, tau = np.broadcast_arrays(*(np.atleast_1d(np.asarray(x, float)) for x in (soil_albedo, rho, tau)))
    bands = soil_albedo.size
    patterns = np.stack([nodes._scatter(1.0, 0.0, 0.0), nodes._scatter(0.0, 1.0, 0.0), nodes._scatter(0.0, 0.0, 1.0)])
    coefficients = np.stack([soil_albedo, (2.0 * rho + tau) / 3.0, (rho + 2.0 * tau) / 3.0], axis=1)
    scatter = np.einsum("bk,kij->bij", coefficients, patterns)
    arrivals = np.broadcast_to(np.asarray(arrivals, float), (bands, 2 * len(nodes.capacity)))
    first = np.einsum("bij,bj->bi", scatter, arrivals)
    if sources is not None:
        first = first + np.broadcast_to(np.asarray(sources, float), first.shape)
    total = np.linalg.solve(np.eye(nodes.size) - scatter @ nodes.arrive.T, first[..., None])[..., 0]
    absorptivity = np.where(nodes.soil, 1.0 - soil_albedo[:, None], 1.0 - rho[:, None] - tau[:, None])

    def per_node(incident):
        return (incident[:, 0::2] + incident[:, 1::2]) * absorptivity
    return per_node(arrivals), per_node(total @ nodes.arrive)


def shortwave_multiple(kind, lai, shape, rho, tau, soil_rho, vza, vaa, sza, saa, row_azimuth=0.0,
                       direct=1.0, diffuse=0.0):
    """Beam orders >= 2 plus all orders of an isotropic sky, with bi-Lambertian leaf scattering."""
    if lai <= 0.0:
        return float(soil_rho) * float(diffuse)
    total = 0.0
    if diffuse > 0.0:
        total += diffuse * float(sky_response(kind, lai, shape, rho, tau, soil_rho, vza, vaa, row_azimuth)[0].sum())
    if direct > 0.0 and sza < 90.0:
        total += direct * float(beam_response(kind, lai, shape, rho, tau, soil_rho, vza, vaa, sza, saa,
                                              row_azimuth)[0].sum())
    return total
