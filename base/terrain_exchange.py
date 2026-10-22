"""Terrain exchange geometry restored from GitHub tirt commit 86c3829."""
from functools import lru_cache

import numpy as np


def _build_hemisphere(n_zenith=None):
    """Zenith nodes (as tan) and sin(theta) d(theta) weights of the upper hemisphere.

    The integrands depend on zenith only, so no azimuth nodes are needed (all uses are normalised).
    None keeps the original 0..89 degree grid of 90 nodes; an integer uses that many midpoint nodes.
    """
    rd = np.pi / 180.0
    if n_zenith is None:
        hza = np.linspace(0, 89, 90)
        dangle = np.pi / 2.0 / 90
    else:
        hza = (np.arange(n_zenith) + 0.5) * 90.0 / n_zenith
        dangle = np.pi / 2.0 / n_zenith
    tantemp = np.tan(hza * rd)
    fweiplus = np.sin(hza * rd) * dangle
    fweiplusSum = np.sum(1.0 * fweiplus)
    tantemp.flags.writeable = False
    fweiplus.flags.writeable = False
    return tantemp, fweiplus, fweiplusSum


_HEMISPHERE = _build_hemisphere(9)


def _mount_projection(dh, radius2, density2, tantemp, outside=False):
    # outside=True: the origin lies outside the cone's cross-section, so only the part beyond it occludes
    if dh <= 0:
        return np.zeros_like(tantemp)
    L2_v = dh * tantemp
    ind = L2_v < radius2
    L2_v[ind] = radius2
    gamma2_v = np.arcsin(radius2 / L2_v)
    offset = -np.pi / 2.0 if outside else np.pi / 2.0
    projtemp_mount = density2 * (1.0/np.tan(gamma2_v) + gamma2_v + offset) * radius2 * radius2
    if outside:
        projtemp_mount[L2_v <= radius2] = 0
    return projtemp_mount


def _hemisphere_sums(shapes, n_part):
    shapes = np.asarray(shapes)
    return _hemisphere_sums_cached(shapes.tobytes(), shapes.shape, shapes.dtype.str, n_part)


def _tilted_view_weights(beta, theta):
    """Lambertian view-factor weight of zenith theta (above: up, below: down) for a facet of slope beta.

    Integrates max(n.omega, 0) over azimuth: a = sin(beta) sin(theta), b = +-cos(beta) cos(theta).
    """
    a = np.sin(beta) * np.sin(theta)
    out = []
    for b in (np.cos(beta) * np.cos(theta), -np.cos(beta) * np.cos(theta)):
        phi0 = np.arccos(np.clip(-b / np.maximum(a, 1e-300), -1.0, 1.0))
        azimuthal = np.where(b >= a, 2 * np.pi * b, np.where(b <= -a, 0.0, 2 * (b * phi0 + a * np.sin(phi0))))
        w = np.sin(theta) * azimuthal
        out.append(w / np.sum(w))
    return out


_FACET_PHI = np.linspace(0.0, 2 * np.pi, 10, endpoint=False)
_SHARE_THETA, _SHARE_PHI = np.meshgrid((np.arange(12) + 0.5) * np.pi / 12, (np.arange(24) + 0.5) * np.pi / 12, indexing="ij")


def _unit(zenith, azimuth):
    return np.stack([np.sin(zenith) * np.cos(azimuth), np.sin(zenith) * np.sin(azimuth), np.cos(zenith) * np.ones_like(azimuth)], -1)


def _cone_normals(beta):
    return _unit(np.full(_FACET_PHI.size, beta), _FACET_PHI)


def _self_lit_fraction(beta, v, s):
    """Exact fraction of the visible (projected) area of a cone side with slope beta whose facets face the sun."""
    normals = _cone_normals(beta)
    w = np.maximum(normals @ v, 0.0)
    return np.sum(w * ((normals @ s) > 0)) / np.sum(w)


def facet_emitter_sunlit_share(emit_normals, emit_area, s, recv_normals, recv_weight, hit_weight):
    """Self-lit share of the facets that irradiate a receiver population.

    A receiver looking along omega (weight max(n_r.omega,0) * P(hit facet)) sees emitter facets facing back
    (weight area * max(-n_e.omega,0)); only those with n_e.s > 0 are sunlit. Validated against Monte Carlo.
    """
    omega = _unit(_SHARE_THETA, _SHARE_PHI).reshape(-1, 3)
    domega = (np.sin(_SHARE_THETA) * hit_weight(_SHARE_THETA)).reshape(-1)
    seen = np.maximum(-(omega @ emit_normals.T), 0.0) * emit_area
    lit = seen @ ((emit_normals @ s) > 0).astype(float) / np.maximum(seen.sum(1), 1e-300)
    weight = (np.maximum(omega @ recv_normals.T, 0.0) @ recv_weight) * domega
    return float(np.sum(weight * lit) / np.sum(weight))


@lru_cache(maxsize=32)
def _terrain_share_geometry(raw, shape, dtype):
    """Sun- and view-independent parts of the facet sunlit-share integrals (cached per terrain shapes)."""
    shapes = np.frombuffer(raw, dtype=dtype).reshape(shape)
    poccupied = np.sum(shapes[:, 2] * np.pi * shapes[:, 1] ** 2)
    normals = np.concatenate([_cone_normals(np.arctan(h / r)) for h, r, d in shapes])
    area = np.concatenate([np.full(_FACET_PHI.size, d * np.pi * r * np.hypot(h, r)) for h, r, d in shapes])
    omega = _unit(_SHARE_THETA, _SHARE_PHI).reshape(-1, 3)
    theta = _SHARE_THETA.reshape(-1)
    tan_theta = np.abs(np.tan(theta))
    upper = theta < np.pi / 2

    seen = np.maximum(-(omega @ normals.T), 0.0) * area
    seen_norm = seen / np.maximum(seen.sum(1, keepdims=True), 1e-300)
    recv_view = np.maximum(omega @ normals.T, 0.0)

    ground_projection = sum(_mount_projection(h, r, d, tan_theta, outside=True) for h, r, d in shapes)
    ground_hit = np.where(upper, 1 - np.exp(-ground_projection / (1 - poccupied)), 0.0)
    ground_weight = np.maximum(omega[:, 2], 0.0) * np.sin(theta) * ground_hit

    facet_weights = []
    for h, r, d in shapes:
        # visible facet points sit at mean height h/3 (band weight 1 - z/h)
        z = h / 3.0
        up = sum(_mount_projection(h2 - z, r2 * (h2 - z) / h2, d2, tan_theta, outside=True) for h2, r2, d2 in shapes if h2 > z)
        down = sum(_mount_projection(min(h2, z), r2, d2, tan_theta) for h2, r2, d2 in shapes)
        facet_weights.append(np.sin(theta) * np.where(upper, 1 - np.exp(-up), 1 - np.exp(-down)))
    return normals, area, seen_norm, recv_view, ground_weight, tuple(facet_weights)


def _terrain_sunlit_shares(shapes, vza, vaa, sza, saa):
    """Sunlit share of facet emission reaching the sensor via facets (s_ww) and via the ground (s_ws).

    Same integrals as facet_emitter_sunlit_share, with the shape-only terms cached.
    """
    rd = np.pi / 180.0
    v = _unit(np.asarray(vza * rd), np.asarray(vaa * rd))
    s = _unit(np.asarray(sza * rd), np.asarray(saa * rd))
    shapes = np.ascontiguousarray(shapes, dtype=float)
    normals, area, seen_norm, recv_view, ground_weight, facet_weights = _terrain_share_geometry(shapes.tobytes(), shapes.shape, shapes.dtype.str)
    lit = seen_norm @ ((normals @ s) > 0).astype(float)

    s_ww, w_sum = 0.0, 0.0
    for k in range(len(shapes)):
        sel = slice(k * _FACET_PHI.size, (k + 1) * _FACET_PHI.size)
        recv = area[sel] * np.maximum(normals[sel] @ v, 0.0)
        if recv.sum() <= 0:
            continue
        weight = (recv_view[:, sel] @ recv) * facet_weights[k]
        s_ww += recv.sum() * float(np.sum(weight * lit) / np.sum(weight))
        w_sum += recv.sum()
    s_ww = s_ww / w_sum
    s_ws = float(np.sum(ground_weight * lit) / np.sum(ground_weight))
    # mutual shadowing of facet points: band-averaged exp(-lambda*A_ext*(1-z/h)^2) = (1-exp(-a))/a
    a = float(sum(_mount_projection(h, r, d, np.atleast_1d(np.tan(sza * rd)), outside=True)[0] for h, r, d in shapes))
    mutual = (1 - np.exp(-a)) / a if a > 1e-12 else 1.0
    return s_ww * mutual, s_ws * mutual


@lru_cache(maxsize=64)
def _hemisphere_sums_cached(raw, shape, dtype, n_part):
    """View-independent upward/downward hemisphere sums of the scattering term."""
    shapes = np.frombuffer(raw, dtype=dtype).reshape(shape)
    tantemp, fweiplus, fweiplusSum = _HEMISPHERE
    up_sums, down_sums = [], []
    theta = np.arctan(tantemp)
    for shape1 in shapes:
        height1 = shape1[0]
        dheight = height1 / n_part
        w_up, w_down = _tilted_view_weights(np.arctan(height1 / shape1[1]), theta)
        up_row, down_row = [], []
        for kh in range(n_part):
            heighttemp = dheight * (kh + 0.5)
            projStemp = 0
            for shape2 in shapes:
                radius2, height2, density2 = shape2[1], shape2[0], shape2[2]
                if (height2 - heighttemp > 0):
                    dh = height2 - heighttemp
                    projStemp = projStemp + _mount_projection(dh, radius2 * dh / height2, density2, tantemp, outside=True)
            up_row.append(np.sum(1.0 * np.exp(-(projStemp)) * w_up))
            projStemp = 0
            for shape2 in shapes:
                radius2, height2, density2 = shape2[1], shape2[0], shape2[2]
                if (height2 - heighttemp < 0):
                    projStemp = projStemp + _mount_projection(height2, radius2, density2, tantemp)
                if (height2 - heighttemp > 0):
                    projStemp = projStemp + _mount_projection(heighttemp, radius2, density2, tantemp)
            down_row.append(np.sum(1.0 * np.exp(-(projStemp)) * w_down))
        up_sums.append(tuple(up_row))
        down_sums.append(tuple(down_row))
    return tuple(up_sums), tuple(down_sums)

###----------------------------------
### ���￼����ɽ�ر������໥�ڵ�
###----------------------------------


def projection_overlap(tv, ts, up):
    """Overlap ratio of view/sun shadow projections; 0 when both directions are vertical."""
    tv = np.asarray(tv, dtype=float)
    ts = np.asarray(ts, dtype=float)
    den = tv + ts
    num = np.sqrt(np.maximum(tv * tv + ts * ts - 2 * tv * ts * up, 0.0))
    result = np.divide(num, den, out=np.zeros(np.broadcast(num, den).shape), where=den > 0)
    return float(result) if result.ndim == 0 else result

class TerrainExchange:
    def scatter_terms(self, kangle):
        """Geometric single-scattering terms shared by Terrain and Terrain_Veg."""
        shapes = self.shapes
        Es = self.Es
        n_shape = self.n_shape
        n_dem = self.n_dim
        n_part = self.n_part
        vza = self.vza_[kangle]
        vaa = self.vaa_[kangle]
        sza = self.sza
        saa = self.saa
        raa = np.abs(vaa - saa) % 360.0

        if raa > 180:raa = 360-raa


        rd = np.pi / 180.0
        ui = np.cos(sza * rd)
        uv = np.cos(vza * rd)
        si = np.sin(sza * rd)
        sv = np.sin(vza * rd)
        up = np.cos(raa * rd)
        tantv = np.tan(vza * rd)
        tants = np.tan(sza * rd)

        ####----------------------------------------------
        ### PLANE surface
        ###-----------------------------------------------
        height1r = 0
        projv = 0
        projs = 0
        poccupied = 0
        projvs = 0
        for kshape2 in range(n_shape):
            shape2 = shapes[kshape2]
            radius2 = shape2[1]
            height2 = shape2[0]
            density2 = shape2[2]
            poccupied = poccupied + density2 * np.pi * (radius2*radius2)
            dh = height2 - height1r
            alpha2 = np.arctan(radius2/height2)

            L2_v = dh * tantv
            if L2_v < radius2: L2_v = radius2*1.0
            theta2_v = np.arctan(L2_v*np.tan(alpha2)/radius2)
            gamma2_v = np.arcsin(radius2/L2_v)
            if dh > 0 and L2_v > radius2:
                projv_mount = density2 * (1.0/np.tan(gamma2_v) + gamma2_v - np.pi/2.0)*radius2*radius2
                projv = projv + projv_mount


            L2_s = dh * tants
            if L2_s < radius2: L2_s = radius2*1.0
            theta2_s = np.arctan(L2_s*np.tan(alpha2)/radius2)
            gamma2_s = np.arcsin(radius2/L2_s)
            if dh > 0 and L2_s > radius2:
                projs_mount = density2 * (1.0/np.tan(gamma2_s) + gamma2_s - np.pi/2.0)*radius2*radius2
                projs = projs + projs_mount

        # ground points lie outside the cone bases: extension-only projections, as in the direct part
        projv = projv/(1-poccupied)
        projs = projs/(1-poccupied)
        Overlapping = projection_overlap(tantv, tants, up)
        projvs = projv + projs * Overlapping
        pPlaneV = np.exp(-projv)*(1-poccupied)
        pPlaneS = np.exp(-projs)*(1-poccupied)
        pPlaneV_sunlit = np.exp(-projvs)*(1-poccupied)
        pPlaneV_shaded = pPlaneV - pPlaneV_sunlit

        if2 = 1.0
        tantemp, fweiplus, fweiplusSum = _HEMISPHERE
        up_sums, down_sums = _hemisphere_sums(shapes, n_part)

        alpha2 = np.arctan(shapes[:, 1] / shapes[:, 0])

        L2_v = shapes[:,0] * tantv
        ind = L2_v < shapes[:, 1]
        L2_v[ind] = shapes[ind, 1]
        theta2_v = np.arctan(L2_v * np.tan(alpha2) / shapes[:, 1])
        gamma2_v = np.arcsin(shapes[:, 1] / L2_v)
        projv_mount = (shapes[:, 2] * (1.0/np.tan(gamma2_v) + gamma2_v + np.pi / 2.0) * shapes[:, 1] * shapes[:, 1])
        tempSv = np.sum(projv_mount)
        f0 = np.exp(-tempSv)
        i0v = 1-f0


        eu = 0.0
        ed = 0.0
        log_weights = []
        up_terms = []
        down_terms = []

        for kshape1 in range(n_shape):
            shape1 = shapes[kshape1]

            radius1 = shape1[1]
            height1 = shape1[0]
            density1 = shape1[2]
            height1r = height1
            dheight = height1 / n_part
            projv = 0
            projs = 0

            alpha1 = np.arctan(radius1 / height1)
            L1_v = height1 * tantv
            if L1_v < radius1: L1_v = radius1 * 1.0
            theta1_v = np.arctan(L1_v * np.tan(alpha1) / radius1)
            gamma1_v = np.arcsin(radius1 / L1_v)
            proj1 = (1.0/np.tan(gamma1_v) + gamma1_v + np.pi / 2.0) * radius1 * radius1
            # up/down split of a facet with slope 90-alpha1: unobstructed view factors (1 +- sin(alpha1))/2
            au = 0.5 * (1 + np.sin(alpha1))
            ad = 0.5 * (1 - np.sin(alpha1))

            for kh in range(n_part):
                heighttemp = dheight * (kh + 0.5)
                # first-hit band of the cone above heighttemp: projection scales with (1 - z/h)^2
                tempdS = density1 * dheight * 2 * (1 - heighttemp / height1) / height1 * proj1

                projStemp = 0
                proj_roof = 0
                tempSv = 0
                for kshape2 in range(n_shape):
                    shape2 = shapes[kshape2]
                    radius2 = shape2[1]
                    height2 = shape2[0]
                    density2 = shape2[2]

                    if (height2 - heighttemp > 0):

                        dh = height2 - heighttemp
                        radius2t = radius2 * dh / height2
                        L2_v = dh * tantv
                        if L2_v < radius2t: L2_v = radius2t * 1.0
                        gamma2_v = np.arcsin(radius2t / L2_v)
                        projv_mount = density2 * (1.0/np.tan(gamma2_v) + gamma2_v + np.pi / 2.0) * radius2t * radius2t
                        tempSv = tempSv + projv_mount

                log_weights.append(np.log(tempdS) - tempSv)
                up_terms.append(au * up_sums[kshape1][kh])
                down_terms.append(ad * down_sums[kshape1][kh])

        # sum(tempdS*exp(-tempSv)) equals i0v analytically; normalising by it keeps the vza->90 limit finite
        weights = np.exp(np.asarray(log_weights) - np.max(log_weights))
        eu = np.sum(weights * np.asarray(up_terms)) / np.sum(weights)
        ed = np.sum(weights * np.asarray(down_terms)) / np.sum(weights)
        if eu + ed > 1.0:
            eu, ed = eu / (eu + ed), ed / (eu + ed)
        p = 1 - eu - ed
        if i0v == 0:
            p = 0
            eu = 0
            ed = 0

        ### ����ֵ���Ȼ��ֵ��ı�����Ȼ��ֵ��뿪���浽����յĸ���
        heighttemp = 0
        ws = 0
        C = 1.05
        projStemp = 0
        tempSv = 0
        for kshape2 in range(n_shape):
            shape2 = shapes[kshape2]
            radius2 = shape2[1]
            height2 = shape2[0]
            density2 = shape2[2]
            if (height2 - heighttemp > 0):
                projStemp = projStemp + _mount_projection(height2 - heighttemp, radius2, density2, tantemp, outside=True)

        ftemp = np.exp(-(projStemp / (1 - poccupied)))
        # a horizontal Lambertian surface weights incoming directions by cos(zenith)
        cosweight = fweiplus * np.cos(np.arctan(tantemp))
        ws = (1-np.sum(if2*(ftemp) * cosweight / np.sum(cosweight))) * pPlaneV

        s_ww, s_ws = _terrain_sunlit_shares(np.asarray(shapes, float), vza, vaa, sza, saa)
        return dict(i0v=i0v, eu=eu, ed=ed, p=p, ws=ws, poccupied=poccupied, pPlaneV=pPlaneV,
                    ground_sunlit=pPlaneS / (1 - poccupied), s_ww=s_ww, s_ws=s_ws)
