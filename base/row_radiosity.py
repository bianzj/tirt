"""Periodic row cross-section radiosity for thermal and shortwave scattering (all orders).

Rows are infinite turbid prisms (width w, gap b, height H, spherical leaves, G = 0.5).
One period of the cross-section is split into soil strips and canopy cells. The angle
out of the cross-section plane is integrated analytically, so with tau the in-plane
optical path a Lambertian soil strip transmits (4/pi) Ki3(tau) and isotropic leaf
emission Ki2(tau) (Bickley functions). The exchange is symmetrised to exact capacity
reciprocity (soil: strip width; canopy cell: 4 k area) and the view sees the nodes
through the same cross-section. Thermal: every node emits with its own area-mean sunlit
fraction. Shortwave: the beam's first scattering leaves the leaves with the spherical
bi-Lambertian phase function, later orders and the sky are isotropic/Lambertian. The
first-order term toward the sensor stays the caller's row geometry.
"""
import functools
import math

import numpy as np

from .row_parametric import G

_PSI, _PSI_W = np.polynomial.legendre.leggauss(64)
_PSI, _PSI_W = 0.25 * np.pi * (_PSI + 1.0), 0.25 * np.pi * _PSI_W
_TAU = np.linspace(0.0, 60.0, 24001)
_DECAY = np.exp(-_TAU[:, None] / np.cos(_PSI))
_KI2 = _DECAY @ (_PSI_W * np.cos(_PSI))
_KI3 = _DECAY @ (_PSI_W * np.cos(_PSI) ** 2) * (4.0 / np.pi)
del _DECAY


def _volume_survival(tau):
    return np.interp(tau, _TAU, _KI2, right=0.0)


def _surface_survival(tau):
    return np.interp(tau, _TAU, _KI3, right=0.0)


def _beam_survival(tau):
    return np.exp(-tau)


_ROW_CELLS, _GAP_CELLS, _LAYERS = 8, 8, 12
_DIRECTIONS = 96
_MAX_PERIODS = 200


class _Cross:
    def __init__(self, lai, width, blank, height):
        self.period, self.height, self.nw, self.nz = width + blank, height, _ROW_CELLS, _LAYERS
        self.k = G * lai * self.period / (width * height)
        gap = np.linspace(width, self.period, _GAP_CELLS + 1)[1:] if blank > 0.0 else np.empty(0)
        self.xe = np.concatenate([np.linspace(0.0, width, _ROW_CELLS + 1), gap])
        self.ze = np.linspace(0.0, height, self.nz + 1)
        self.nx = self.xe.size - 1
        self.n = self.nx + self.nw * self.nz
        dz = height / self.nz
        self.capacity = np.r_[np.diff(self.xe), np.full(self.nw * self.nz, 4.0 * self.k * width / self.nw * dz)]

    def points(self, soil_per_strip, canopy_per_side):
        """Stratified sample points (x, z, owner node): soil strips at z = 0, then canopy cells."""
        u = (np.arange(soil_per_strip) + 0.5) / soil_per_strip
        left, span = self.xe[:-1], np.diff(self.xe)
        xs = (left[:, None] + span[:, None] * u).ravel()
        owner_s = np.repeat(np.arange(self.nx), soil_per_strip)
        v = (np.arange(canopy_per_side) + 0.5) / canopy_per_side
        iz, ix = np.divmod(np.arange(self.nw * self.nz), self.nw)
        xc, zc = np.broadcast_arrays((ix[:, None, None] + v[None, :, None]) * (self.xe[1] - self.xe[0]),
                                     (iz[:, None, None] + v[None, None, :]) * (self.height / self.nz))
        owner_c = np.repeat(self.nx + np.arange(self.nw * self.nz), canopy_per_side ** 2)
        return (xs, np.zeros_like(xs), owner_s), (xc.ravel(), zc.ravel(), owner_c)

    def trace(self, x0, z0, dx, dz):
        """Canopy node and cumulative in-plane optical depth of each segment, plus the exit node.

        (dx, dz) is a unit in-plane direction with dz != 0; the exit node is the soil strip hit
        at z = 0, or -1 when the ray leaves through the canopy top.
        """
        p, count = self.period, x0.size
        t_end = (self.height - z0) / dz if dz > 0.0 else -z0 / dz
        if abs(dx) > 1.0e-12:
            t_end = np.minimum(t_end, _MAX_PERIODS * p / abs(dx))
            n = int(math.ceil(abs(dx) * float(t_end.max()) / p)) + 1
            edges = (self.xe[:-1][None, :] + p * np.arange(-n, n + 1)[:, None]).ravel()
            crossings = [(edges[None, :] - x0[:, None]) / dx]
        else:
            crossings = []
        crossings.append((self.ze[None, :] - z0[:, None]) / dz)
        t = np.concatenate(crossings, axis=1)
        t = np.sort(np.where((t > 1.0e-12) & (t < t_end[:, None]), t, np.inf), axis=1)
        t = np.minimum(t[:, :int(np.isfinite(t).sum(axis=1).max())], t_end[:, None])
        t = np.concatenate([np.zeros((count, 1)), t, t_end[:, None]], axis=1)
        mid = 0.5 * (t[:, :-1] + t[:, 1:])
        ix = np.clip(np.searchsorted(self.xe, np.mod(x0[:, None] + dx * mid, p), side="right") - 1, 0, self.nx - 1)
        iz = np.clip(np.floor((z0[:, None] + dz * mid) / (self.height / self.nz)).astype(int), 0, self.nz - 1)
        inside = ix < self.nw
        node = np.where(inside, self.nx + iz * self.nw + ix, -1)
        depth = np.concatenate([np.zeros((count, 1)), np.cumsum(np.where(inside, self.k * np.diff(t, axis=1), 0.0), axis=1)], axis=1)
        if dz < 0.0:
            landing = np.mod(x0 + dx * t_end, p)
            exit_node = np.clip(np.searchsorted(self.xe, landing, side="right") - 1, 0, self.nx - 1)
        else:
            exit_node = np.full(count, -1)
        return node, depth, exit_node

    def collide(self, x0, z0, dx, dz, survival, scale=1.0):
        """First-collision distribution (points x nodes) of rays along (dx, dz); escapes are dropped."""
        node, depth, exit_node = self.trace(x0, z0, dx, dz)
        alive = survival(depth * scale)
        count, slots = x0.size, self.n + 1
        rows = np.repeat(np.arange(count), node.shape[1])
        index = np.concatenate([rows * slots + np.where(node < 0, self.n, node).ravel(),
                                np.arange(count) * slots + np.where(exit_node < 0, self.n, exit_node)])
        weight = np.concatenate([(alive[:, :-1] - alive[:, 1:]).ravel(), alive[:, -1]])
        return np.bincount(index, weight, count * slots).reshape(count, slots)[:, :self.n]


def _in_plane(zenith, azimuth, row_azimuth):
    """Unit in-plane direction (across-row, up) of a 3-D direction and its in-plane length cos(psi)."""
    theta = math.radians(min(abs(float(zenith)), 89.999))
    across = math.sin(theta) * math.sin(math.radians(float(azimuth) - float(row_azimuth)))
    length = math.hypot(across, math.cos(theta))
    return across / length, math.cos(theta) / length, length


@functools.lru_cache(maxsize=16)
def _exchange(lai, width, blank, height):
    """Reciprocal flux exchange factors F (n x n) of one row period."""
    cross = _Cross(lai, width, blank, height)
    (xs, zs, _), (xc, zc, _) = cross.points(4, 2)
    factors = np.zeros((cross.n, cross.n))
    phi = (np.arange(_DIRECTIONS) + 0.5) * 2.0 * np.pi / _DIRECTIONS
    upward = np.sin(phi) > 0.0
    lambert = np.sin(phi[upward]).sum()
    for angle in phi:
        dx, dz = math.cos(angle), math.sin(angle)
        hits = cross.collide(xc, zc, dx, dz, _volume_survival)
        factors[cross.nx:] += hits.reshape(cross.nw * cross.nz, -1, cross.n).mean(axis=1) / _DIRECTIONS
        if dz > 0.0:
            hits = cross.collide(xs, zs, dx, dz, _surface_survival)
            factors[:cross.nx] += hits.reshape(cross.nx, -1, cross.n).mean(axis=1) * dz / lambert
    exchange = cross.capacity[:, None] * factors
    return cross, 0.5 * (exchange + exchange.T) / cross.capacity[:, None]


@functools.lru_cache(maxsize=64)
def _scattering(lai, width, blank, height, leaf_albedo, soil_albedo):
    """Operator mapping node sources to the light each node scatters onwards (all further orders)."""
    cross, factors = _exchange(lai, width, blank, height)
    albedo = np.r_[np.full(cross.nx, soil_albedo), np.full(cross.n - cross.nx, leaf_albedo)]
    response = np.linalg.solve(np.eye(cross.n) - albedo[:, None] * factors, np.eye(cross.n))
    return cross, albedo, factors, albedo[:, None] * (factors @ response)


def sunlit_fractions(cross, sza, saa, row_azimuth):
    """Area-mean sunlit fraction of every node."""
    if sza >= 90.0:
        return np.zeros(cross.n)
    dx, dz, length = _in_plane(sza, saa, row_azimuth)
    x0, z0, owner = [np.concatenate(part) for part in zip(*cross.points(8, 3))]
    _, depth, _ = cross.trace(x0, z0, dx, dz)
    return np.bincount(owner, np.exp(-depth[:, -1] / length), cross.n) / np.bincount(owner, minlength=cross.n)


def view_weights(cross, vza, vaa, row_azimuth, entries_per_cell=8):
    """Fraction of the view footprint whose first collision is in each node."""
    dx, dz, length = _in_plane(vza, vaa, row_azimuth)
    count = cross.nx * entries_per_cell
    x0 = (np.arange(count) + 0.5) * cross.period / count
    hits = cross.collide(x0, np.full(count, cross.height), -dx, -dz, _beam_survival, 1.0 / length)
    return hits.mean(axis=0)


def _seen(cross, vza, vaa, row_azimuth, soil_view):
    seen = view_weights(cross, vza, vaa, row_azimuth)
    if soil_view is not None:
        soil = np.arange(cross.n) < cross.nx
        for part, target in ((soil, float(soil_view)), (~soil, 1.0 - float(soil_view))):
            total = seen[part].sum()
            seen[part] = seen[part] * target / total if total > 0.0 else 0.0
    return seen


def scatter_emissivity(lai, width, blank, height, emis_leaf, emis_soil, vza, vaa, sza, saa, row_azimuth,
                       soil_view=None):
    """Scattered component weights [soil sunlit, soil shaded, leaf sunlit, leaf shaded] for one direction.

    soil_view rescales the visible soil/leaf nodes to the caller's direct view fractions, so
    direct + scattered stays <= 1 under a black sky.
    """
    if lai <= 0.0:
        return np.zeros(4)
    cross, _, _, operator = _scattering(float(lai), float(width), float(blank), float(height),
                                        1.0 - float(emis_leaf), 1.0 - float(emis_soil))
    sunlit = sunlit_fractions(cross, sza, saa, row_azimuth)
    soil = np.arange(cross.n) < cross.nx
    source = np.zeros((cross.n, 4))
    source[soil, 0], source[soil, 1] = emis_soil * sunlit[soil], emis_soil * (1.0 - sunlit[soil])
    source[~soil, 2], source[~soil, 3] = emis_leaf * sunlit[~soil], emis_leaf * (1.0 - sunlit[~soil])
    return _seen(cross, vza, vaa, row_azimuth, soil_view) @ (operator @ source)


_MU, _MU_W = np.polynomial.legendre.leggauss(16)
_AZIMUTHS = 32


@functools.lru_cache(maxsize=32)
def _first_redistribution(lai, width, blank, height, sza, relative_azimuth):
    """Where light scattered once by sunlit leaves first collides (canopy nodes x nodes).

    Spherical bi-Lambertian leaves: Gamma = rho (sin b - b cos b)/(3 pi) + tau [(sin b - b cos b)/(3 pi)
    + cos b / 3], b the scattering angle. Returns the reflected and transmitted parts, each
    normalised to unit scattered power.
    """
    cross, _ = _exchange(lai, width, blank, height)
    _, (xc, zc, _) = cross.points(4, 2)
    theta, delta = math.radians(sza), math.radians(relative_azimuth)
    incident = -np.array([math.sin(theta) * math.sin(delta), math.sin(theta) * math.cos(delta), math.cos(theta)])
    parts = np.zeros((2, cross.n - cross.nx, cross.n))
    total = np.zeros(2)
    for mu, weight in zip(_MU, _MU_W):
        sine = math.sqrt(1.0 - mu * mu)
        for phi in (np.arange(_AZIMUTHS) + 0.5) * 2.0 * np.pi / _AZIMUTHS:
            out = np.array([sine * math.sin(phi), sine * math.cos(phi), mu])
            cosine = float(np.clip(out @ incident, -1.0, 1.0))
            angle = math.acos(cosine)
            diffuse = (math.sin(angle) - angle * cosine) / (3.0 * np.pi)
            phase = np.array([diffuse, diffuse + cosine / 3.0]) * weight * 2.0 * np.pi / _AZIMUTHS
            total += phase
            length = math.hypot(out[0], out[2])
            hits = cross.collide(xc, zc, out[0] / length, out[2] / length, _beam_survival, 1.0 / length)
            parts += phase[:, None, None] * hits.reshape(cross.nw * cross.nz, -1, cross.n).mean(axis=1)[None]
    return parts / total[:, None, None]


def shortwave_multiple(lai, width, blank, height, leaf_rho, leaf_tau, soil_albedo, vza, vaa, sza, saa,
                       row_azimuth, direct=1.0, diffuse=0.0, soil_view=None):
    """Reflectance factor of the beam's second and higher orders plus all orders of an isotropic sky.

    direct/diffuse are the beam and sky shares of the horizontal irradiance. The beam's first
    scattering leaves the leaves with the spherical bi-Lambertian phase function and the soil
    Lambertian; later orders and the sky are Lambertian/isotropic. The beam's first order
    toward the sensor is the caller's.
    """
    if lai <= 0.0:
        return float(soil_albedo) * float(diffuse)
    leaf_albedo = float(leaf_rho) + float(leaf_tau)
    cross, albedo, factors, operator = _scattering(float(lai), float(width), float(blank), float(height),
                                                   leaf_albedo, float(soil_albedo))
    sky = np.clip(1.0 - factors.sum(axis=1), 0.0, 1.0)
    # The sky term uses node/sky reciprocity, renormalised so the nodes intercept exactly the sky flux.
    sky *= cross.period / float(np.sum(cross.capacity * sky))
    leaving = diffuse * albedo * sky / np.pi
    leaving = leaving + operator @ leaving
    if direct > 0.0 and sza < 90.0:
        beam = direct * view_weights(cross, sza, saa, row_azimuth) * cross.period
        reflected, transmitted = _first_redistribution(float(lai), float(width), float(blank), float(height),
                                                       round(float(sza), 6), round((float(saa) - float(row_azimuth)) % 360.0, 6))
        incident = (float(soil_albedo) * beam[:cross.nx]) @ factors[:cross.nx]
        incident += beam[cross.nx:] @ (float(leaf_rho) * reflected + float(leaf_tau) * transmitted)
        second = albedo * incident / (np.pi * cross.capacity)
        leaving += second + operator @ second
    return float(np.pi * _seen(cross, vza, vaa, row_azimuth, soil_view) @ leaving)
