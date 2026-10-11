"""Compound-terrain component proportions based on TIRT terrain_veg.py."""

import math

import numpy as np

from .physicsF import slope1


def _shade(values):
    values = np.asarray(values, dtype=float).reshape(4)
    return np.asarray([0.0, values[0] + values[1], 0.0, values[2] + values[3]])


def _facet_directions(zenith, azimuth, slope_deg, facet_azimuths):
    """Vectorized solar.slope_direction over facet aspects; slope1 depends on the aspect only via azimuth - aspect."""
    count = facet_azimuths.size
    relative = float(azimuth) - facet_azimuths
    with np.errstate(divide="ignore", invalid="ignore"):
        local_zenith, local_azimuth = slope1(np.full(count, float(zenith)), relative, float(slope_deg), 0.0)
    rd = math.pi / 180.0
    incidence = (
        math.cos(float(slope_deg) * rd) * math.cos(float(zenith) * rd)
        + math.sin(float(slope_deg) * rd) * math.sin(float(zenith) * rd) * np.cos(relative * rd)
    )
    return local_zenith, local_azimuth, incidence >= 0.0


def _shape_values(settings):
    return tuple(
        (max(float(height), 1.0e-6), max(float(radius), 1.0e-6), max(float(density), 0.0))
        for height, radius, density in settings["terrain_shapes"]
    )


def _safe_tangent(angle):
    cosine = math.cos(math.radians(float(angle)))
    if abs(cosine) < 1.0e-8:
        return math.copysign(1.0e8, math.sin(math.radians(float(angle))))
    return math.tan(math.radians(float(angle)))


def _projected_gap(shapes, vza, sza, raa, n_part=10):
    tv = max(_safe_tangent(vza), 0.0)
    ts = max(_safe_tangent(sza), 0.0)
    overlap_denominator = tv + ts
    relative = math.cos(math.radians(raa))
    occupied = sum(density * math.pi * radius * radius for height, radius, density in shapes)
    occupied = min(max(occupied, 0.0), 0.999999)
    plane_v = 0.0
    plane_s = 0.0
    for height, radius, density in shapes:
        for tangent, target in ((tv, "view"), (ts, "sun")):
            length = max(height * tangent, radius)
            gamma = math.asin(min(1.0, radius / length))
            if height > 0.0 and length > radius:
                projected = density * (1.0 / math.tan(gamma) + gamma - math.pi / 2.0) * radius * radius
                if target == "view":
                    plane_v += projected
                else:
                    plane_s += projected
    plane_v /= max(1.0 - occupied, 1.0e-12)
    plane_s /= max(1.0 - occupied, 1.0e-12)
    if overlap_denominator > 1.0e-12:
        overlap = math.sqrt(max(tv * tv + ts * ts - 2.0 * tv * ts * relative, 0.0)) / overlap_denominator
    else:
        overlap = 0.0
    plane_vs = plane_v + plane_s * overlap
    p_plane = math.exp(-min(max(plane_v, -700.0), 700.0)) * (1.0 - occupied)
    p_plane_sun = math.exp(-min(max(plane_vs, -700.0), 700.0)) * (1.0 - occupied)
    return min(max(p_plane, 0.0), 1.0), min(max(p_plane_sun, 0.0), 1.0), occupied


def _mountain_visibility(shapes, vza, sza, raa, n_part=10):
    tv = max(_safe_tangent(vza), 0.0)
    ts = max(_safe_tangent(sza), 0.0)
    overlap_denominator = tv + ts
    relative = math.cos(math.radians(raa))
    p_mount = 0.0
    p_mount_sun = 0.0
    details = []
    for height, radius, density in shapes:
        alpha = math.atan2(radius, height)
        slope = math.degrees(alpha)
        lv = max(height * tv, radius)
        ls = max(height * ts, radius)
        gamma_v = math.asin(min(1.0, radius / lv))
        gamma_s = math.asin(min(1.0, radius / ls))
        projected_v = 0.0
        projected_s = 0.0
        weights = []
        for part in range(n_part):
            z = height * (part + 0.5) / n_part
            dz = height - z
            for tangent, key in ((tv, "view"), (ts, "sun")):
                length = np.maximum(dz * tangent, radius)
                gamma = np.arcsin(np.minimum(1.0, radius / length))
                projected = np.where(
                    (dz > 0.0) & (length > radius),
                    density * (1.0 / np.tan(gamma) + gamma - math.pi / 2.0) * radius * radius,
                    0.0,
                )
                if key == "view":
                    projected_v += float(projected)
                else:
                    projected_s += float(projected)
            radius_at_height = z / height * radius
            ring_weight = (radius_at_height * math.pi * (height / n_part) + radius_at_height * 2.0 + height / n_part)
            weights.append(ring_weight)
        normalizer = max(n_part * sum(weights), 1.0e-12)
        projected_v /= normalizer * max(1.0 - sum(d * math.pi * r * r for _, r, d in shapes), 1.0e-12)
        projected_s /= normalizer * max(1.0 - sum(d * math.pi * r * r for _, r, d in shapes), 1.0e-12)
        if overlap_denominator > 1.0e-12:
            overlap = math.sqrt(max(tv * tv + ts * ts - 2.0 * tv * ts * relative, 0.0)) / overlap_denominator
        else:
            overlap = 0.0
        projected_vs = projected_v + projected_s * overlap
        gap_v = math.exp(-min(max(projected_v, -700.0), 700.0))
        gap_vs = math.exp(-min(max(projected_vs, -700.0), 700.0))
        if lv <= radius:
            visible_area = density * math.pi * radius * radius * gap_v
            sun_area = density * math.pi * radius * radius * gap_vs
            if sza <= slope:
                sun_area = density * math.pi * radius * radius * gap_v
        else:
            area = (1.0 / math.tan(gamma_v) + gamma_v + math.pi / 2.0) * radius * radius
            visible_area = density * gap_v * area
            sun_area = density * gap_vs * area
            if sza <= slope:
                sun_area = density * gap_v * area
        p_mount += visible_area
        p_mount_sun += sun_area
        details.append((alpha, visible_area, sun_area))
    return p_mount, p_mount_sun, details


def compound_terrain_proportions(settings, sza, saa, vza, vaa, vegetation_proportions, sunlit=True):
    """Return soil/vegetation sunlit/shaded fractions for TIRT-like terrain."""
    shapes = _shape_values(settings)
    n_part = 10
    sza = min(max(float(sza), 0.0), 89.9)
    vza = min(max(float(vza), 0.0), 89.9)
    saa = float(saa)
    vaa = float(vaa)
    raa = abs(vaa - saa) % 360.0
    if raa > 180.0:
        raa = 360.0 - raa
    p_plane, p_plane_sun, occupied = _projected_gap(shapes, vza, sza, raa, n_part)
    p_mount_raw, p_mount_sun_raw, mountain_details = _mountain_visibility(shapes, vza, sza, raa, n_part)
    p_mount = max(1.0 - p_plane, 0.0)
    scale = p_mount / max(p_mount_raw, 1.0e-12)
    p_mount_sun = min(max(p_mount_sun_raw * scale, 0.0), p_mount)
    p_mount_shade = p_mount - p_mount_sun

    plane = np.asarray(vegetation_proportions(sza, saa, vza, vaa, True), dtype=float).reshape(4)
    if not sunlit:
        plane = _shade(plane)
    if not sunlit:
        p_plane_sun = 0.0
    plane_sun = np.asarray([p_plane_sun * plane[0], 0.0, p_plane_sun * plane[2], 0.0])
    plane_shade = np.asarray([
        0.0,
        p_plane_sun * plane[1] + (p_plane - p_plane_sun) * (plane[0] + plane[1]),
        0.0,
        p_plane_sun * plane[3] + (p_plane - p_plane_sun) * (plane[2] + plane[3]),
    ])

    mountain = np.zeros(4, dtype=float)
    mountain_sun = np.zeros(4, dtype=float)
    for alpha, visible_area, sun_area in mountain_details:
        if visible_area <= 0.0:
            continue
        azimuths = np.arange(0.0, 360.0, 5.0)
        cos_view = (
            math.cos(alpha) * math.cos(math.radians(vza))
            + math.sin(alpha) * math.sin(math.radians(vza))
            * np.cos(np.radians(vaa - azimuths))
        )
        visible = cos_view >= 0.0
        if not np.any(visible):
            continue
        weights = np.maximum(cos_view[visible], 0.0)
        weights /= max(float(weights.sum()), 1.0e-12)
        facet_azimuths = azimuths[visible]
        view_z, view_a, _ = _facet_directions(vza, vaa, math.degrees(alpha), facet_azimuths)
        sun_z, sun_a, sun_lit = _facet_directions(sza, saa, math.degrees(alpha), facet_azimuths)
        facet_visible_total = np.zeros(4, dtype=float)
        facet_sun_total = np.zeros(4, dtype=float)
        for index, weight in enumerate(weights):
            facet_sunlit = bool(sun_lit[index])
            facet = np.asarray(
                vegetation_proportions(
                    float(sun_z[index]), float(sun_a[index]), float(view_z[index]), float(view_a[index]),
                    facet_sunlit and sunlit,
                ),
                dtype=float,
            ).reshape(4)
            facet_visible_total += weight * facet
            if facet_sunlit and sunlit:
                facet_sun_total += weight * facet
        mountain += visible_area * facet_visible_total
        if sun_area > 0.0 and sunlit:
            mountain_sun += sun_area * facet_sun_total
    mountain *= scale
    mountain_sun = np.minimum(np.maximum(mountain_sun * scale, 0.0), mountain)
    mountain_shade = mountain - mountain_sun

    result = plane_sun + plane_shade + mountain_sun + mountain_shade
    result = np.maximum(result, 0.0)
    total = float(result.sum())
    return result / total if total > 1.0e-12 else np.asarray([1.0, 0.0, 0.0, 0.0])
