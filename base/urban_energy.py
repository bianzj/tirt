"""Conservative city energy exchange with a canopy/soil courtyard boundary.

The five macro nodes are courtyard plus four facade bearings. The courtyard's
reflection/emission and component absorption use the same canopy node operator;
canopy leaves add their actual area rather than replacing visible street area.
"""
import math

import numpy as np

from . import node_radiosity as nr
from . import tirteb_radiation as tr
from .optical_urban import building_shape, obstruction, wall_exchange

SIGMA = 5.67e-8


def geometry(settings, state):
    shape = building_shape(settings)
    length, width, height, density, bearing = shape
    roof, ground = density * length * width, 1. - density * length * width
    walls = density * height * np.array([length, width, length, width])
    if min(walls) > 1.e-14:
        factors, sky, walls = wall_exchange(shape, 0.)
    else:
        factors, sky = np.zeros((5, 5)), np.ones(5)
    sza, saa = float(state.get("sza", 90.)), float(state.get("saa", 0.))
    lit = bool(state.get("sunlit_surface", True)) and sza < 90.
    incidence = np.maximum(np.cos(np.deg2rad(saa - bearing - np.arange(4) * 90.)), 0.)
    illuminated = incidence > 1.e-12 if lit and sza > 1.e-10 else np.zeros(4, bool)
    building_sun = math.exp(-obstruction(shape, min(sza, 89.9), saa)) if lit else 0.
    projected = float(walls @ incidence)
    # Allocate only the direct power intercepted above the courtyard. The roof
    # and courtyard+walls then intercept exactly one horizontal beam unit.
    wall_beam = ground * (1. - building_sun) * incidence / projected if projected > 0. else np.zeros(4)
    return dict(shape=shape, roof=roof, ground=ground, areas=np.r_[ground, walls], walls=walls,
                factors=factors, sky=sky, sun=illuminated, beam=wall_beam,
                building_sun=building_sun, lit=lit, sza=sza, saa=saa,
                roof_sunlit_fraction=float(lit),
                wall_sunlit_fraction=float(walls @ illuminated / walls.sum()) if walls.sum() else 0.,
                street_building_sunlit_fraction=building_sun)


def _shares(result, geo):
    for key in ("roof_sunlit_fraction", "wall_sunlit_fraction", "street_building_sunlit_fraction"):
        result[key] = geo[key]
    result.update(urban_street_svf=float(geo["sky"][0]),
                  urban_wall_svf=float(np.average(geo["sky"][1:], weights=geo["walls"])) if geo["walls"].sum() else 1.,
                  urban_leaf_svf=float(geo["sky"][0]))


def _facade_means(values, geo):
    """Area means of the illuminated/unilluminated facade groups."""
    values = np.atleast_2d(values)
    result = []
    for mask in (geo["sun"], ~geo["sun"]):
        weights = geo["walls"] * mask
        if weights.sum() <= 1.e-14:
            weights = geo["walls"]
        # Zero-area diagnostic slots use the finite open-boundary limit. They
        # carry no scene power and must not introduce a spurious cold wall.
        result.append(values @ (weights / weights.sum()) if weights.sum() else values.mean(axis=1))
    return np.stack(result, axis=-1)


def _canopy(settings, geo):
    layout = tr._scene(settings, geo["sza"], geo["saa"], geo["lit"])
    if layout is None:
        return None
    nodes, beam, soil_area, leaf_area, sun_area = layout
    return nodes, beam, soil_area, leaf_area, sun_area * geo["building_sun"]


def _component_means(power, canopy):
    nodes, _, soil_area, leaf_area, sun_area = canopy
    soil = tr._means(np.atleast_2d(power), soil_area, sun_area, nodes.soil)
    leaf = tr._means(np.atleast_2d(power), leaf_area, sun_area, ~nodes.soil)
    return np.concatenate((soil, leaf), axis=-1)


def _thermal_boundary(settings, state, geo):
    es = float(settings.get("urban_emissivity_street", .955))
    el = float(settings.get("leaf_emissivity", .985))
    canopy = _canopy(settings, geo)
    temps = np.array([float(state.get(k, 293.15)) ** 4 for k in ("Tss", "Tsh", "Tls", "Tlh")])
    if canopy is None:
        share = geo["building_sun"]
        emitted = es * SIGMA * (share * temps[0] + (1. - share) * temps[1])
        return 1. - es, emitted, np.array([es, es, 0., 0.]), np.zeros(4)
    nodes, _, soil_area, leaf_area, sun_area = canopy
    share = np.divide(sun_area, soil_area + leaf_area, out=np.zeros_like(sun_area), where=soil_area + leaf_area > 0.)
    emit = np.where(nodes.soil, es * (share * temps[0] + (1. - share) * temps[1]),
                    el * (share * temps[2] + (1. - share) * temps[3])) * SIGMA
    sources = nodes.emission(emit[:, None] / math.pi)[:, 0]
    first, scattered = nr.absorbed(nodes, 1. - es, 1. - el, 0., nr.sky_arrivals(nodes))
    sky_abs = first + scattered
    first, scattered = nr.absorbed(nodes, 1. - es, 1. - el, 0., np.zeros(2 * len(nodes.capacity)), sources)
    internal = first + scattered
    # Kirchhoff-consistent effective boundary: emitted power escaping the local
    # canopy, and reflection of a unit isotropic field over that same boundary.
    reflectance = 1. - float(sky_abs.sum())
    emitted = float(nodes.capacity @ emit - internal.sum())
    return reflectance, emitted, _component_means(sky_abs, canopy)[0], _component_means(internal, canopy)[0]


def longwave_inputs(settings, state, sky):
    geo = geometry(settings, state)
    sky = max(float(sky), 0.)
    ground_rho, ground_emit, response, internal = _thermal_boundary(settings, state, geo)
    ew = float(settings.get("urban_emissivity_wall", .93))
    er = float(settings.get("urban_emissivity_roof", .95))
    wall_t4 = np.where(geo["sun"], float(state.get("Twall_sunlit", 293.15)) ** 4,
                      float(state.get("Twall_shaded", 293.15)) ** 4)
    rho = np.r_[ground_rho, np.full(4, 1. - ew)]
    emitted = np.r_[ground_emit, ew * SIGMA * wall_t4]
    radiosity = np.linalg.solve(np.eye(5) - rho[:, None] * geo["factors"], emitted + rho * geo["sky"] * sky)
    incident = geo["factors"] @ radiosity + geo["sky"] * sky
    courtyard = internal + response * incident[0]
    wall_means = _facade_means(ew * incident[1:], geo)[0]
    result = {}
    for index, suffix in enumerate(("sunlit", "shaded")):
        result[f"roof_{suffix}"] = er * sky
        result[f"wall_{suffix}"] = wall_means[index]
        result[f"street_{suffix}"] = float(courtyard[index])
        result[f"leaf_{suffix}"] = float(courtyard[index + 2])
        result[f"LW_down_roof_{suffix}"] = sky
        result[f"LW_down_wall_{suffix}"] = float(_facade_means(incident[1:], geo)[0, index])
        result[f"LW_down_street_{suffix}"] = float(incident[0])
        result[f"LW_down_leaf_{suffix}"] = float(incident[0])
    roof_share = geo["roof_sunlit_fraction"]
    roof_t4 = roof_share * float(state.get("Troof_sunlit", 293.15)) ** 4 + (1. - roof_share) * float(state.get("Troof_shaded", 293.15)) ** 4
    roof_out = er * SIGMA * roof_t4 + (1. - er) * sky
    result["urban_LW_escaped"] = float(geo["roof"] * roof_out + geo["areas"] @ (geo["sky"] * radiosity))
    result["urban_LW_net_ground"] = sky - result["urban_LW_escaped"]
    _shares(result, geo)
    return result


def _spectrum(settings, key, fallback, bands):
    values = np.asarray(settings.get(key, ()), float)
    if values.size >= bands:
        return np.clip(values[:bands], 0., 1.)
    return np.full(bands, np.clip(float(fallback), 0., 1.))


def shortwave_inputs(settings, state, netrad):
    """Absorbed city SW per component area, plus escaped power per ground area."""
    geo, result = geometry(settings, state), {}
    beam = np.asarray(netrad.get("directVrad_incident_spectrum", ()), float)
    diffuse = np.asarray(netrad.get("diffuseVrad_incident_spectrum", ()), float)
    if beam.size != tr.BANDS or diffuse.size != tr.BANDS:
        beam = np.array([float(netrad.get("directVrad_incident", 0.))])
        diffuse = np.array([float(netrad.get("diffuseVrad_incident", netrad.get("diffuseVrad_soil", 0.)))])
    beam = np.maximum(beam, 0.) * float(geo["lit"])
    diffuse = np.maximum(diffuse, 0.)
    bands = len(beam)
    rr, rw, rs = (_spectrum(settings, f"urban_reflectance_{component}_spectrum",
                            1. - float(settings.get(f"urban_absorption_{component}", default)), bands)
                  for component, default in (("roof", .95), ("wall", .92), ("street", .955)))
    lr = _spectrum(settings, "leaf_reflectance_spectrum", settings.get("leaf_reflectance", .015), bands)
    lt = np.minimum(_spectrum(settings, "leaf_transmittance_spectrum", settings.get("leaf_transmittance", 0.), bands), 1. - lr)
    canopy = _canopy(settings, geo)
    if canopy is None:
        sky_abs, first_abs, scatter_abs = (1. - rs)[:, None], (1. - rs)[:, None], np.zeros((bands, 1))
    else:
        nodes, first_beam, _, _, _ = canopy
        a, b = nr.absorbed(nodes, rs, lr, lt, nr.sky_arrivals(nodes))
        sky_abs = a + b
        first_abs, scatter_abs = nr.absorbed(nodes, rs, lr, lt, nr.from_above(nodes, first_beam))
    rho_ground = 1. - sky_abs.sum(axis=1)
    beam_escape = 1. - (first_abs + scatter_abs).sum(axis=1)
    rho = np.column_stack((rho_ground, np.repeat(rw[:, None], 4, axis=1)))
    first_sources = np.column_stack((beam_escape * beam * geo["building_sun"],
                                    rw[:, None] * beam[:, None] * geo["beam"]))
    sources = first_sources + rho * diffuse[:, None] * geo["sky"]
    radiosity = np.linalg.solve(np.eye(5) - rho[:, :, None] * geo["factors"], sources[..., None])[..., 0]
    incident_diffuse = radiosity @ geo["factors"].T + diffuse[:, None] * geo["sky"]
    direct_node = beam[:, None] * geo["building_sun"] * first_abs
    diffuse_node = beam[:, None] * geo["building_sun"] * scatter_abs + incident_diffuse[:, :1] * sky_abs
    if canopy is None:
        direct = np.column_stack((beam * (1. - rs), np.zeros((bands, 3))))
        courtyard_diffuse = np.column_stack((diffuse_node[:, 0], diffuse_node[:, 0], np.zeros((bands, 2))))
    else:
        nodes, _, _, _, sun_area = canopy
        direct = np.zeros((bands, 4))
        for column, mask in ((0, nodes.soil), (2, ~nodes.soil)):
            if sun_area[mask].sum() > 1.e-14:
                direct[:, column] = direct_node[:, mask].sum(axis=1) / sun_area[mask].sum()
        courtyard_diffuse = _component_means(diffuse_node, canopy)
    wall_direct = _facade_means((1. - rw[:, None]) * beam[:, None] * geo["beam"], geo)
    wall_diffuse = _facade_means((1. - rw[:, None]) * incident_diffuse[:, 1:], geo)
    for index, suffix in enumerate(("sunlit", "shaded")):
        result[f"SW_direct_roof_{suffix}"] = float(((1. - rr) * beam).sum()) if index == 0 else 0.
        result[f"SW_diffuse_roof_{suffix}"] = float(((1. - rr) * diffuse).sum())
        result[f"SW_direct_wall_{suffix}"] = float(wall_direct[:, index].sum()) if index == 0 else 0.
        result[f"SW_diffuse_wall_{suffix}"] = float(wall_diffuse[:, index].sum())
        for component, column in (("street", index), ("leaf", index + 2)):
            result[f"SW_direct_{component}_{suffix}"] = float(direct[:, column].sum())
            result[f"SW_diffuse_{component}_{suffix}"] = float(courtyard_diffuse[:, column].sum())
    escaped = geo["roof"] * float((rr * (beam + diffuse)).sum()) + float((radiosity * geo["sky"] * geo["areas"]).sum())
    result.update(urban_SW_escaped=escaped, urban_SW_absorbed=float((beam + diffuse).sum()) - escaped,
                  urban_SW_incoming=float((beam + diffuse).sum()),
                  shortwave_albedo=escaped / max(float((beam + diffuse).sum()), 1.e-12))
    for component, offset in (("soil", 0), ("leaf", 2)):
        result[f"directVrad_{component}"] = float(direct[:, offset].sum())
        for index, suffix in enumerate(("sunlit", "shaded")):
            result[f"diffuseVrad_{component}_{suffix}"] = float(courtyard_diffuse[:, offset + index].sum())
        if canopy is None:
            result[f"diffuseVrad_{component}"] = float(diffuse_node.sum()) if component == "soil" else 0.
        else:
            nodes, _, soil_area, leaf_area, _ = canopy
            mask = nodes.soil if component == "soil" else ~nodes.soil
            area = soil_area.sum() if component == "soil" else leaf_area.sum()
            result[f"diffuseVrad_{component}"] = float(diffuse_node[:, mask].sum() / area) if area else 0.
    if bands == tr.BANDS:
        wavelengths = np.asarray(settings.get("wavelengths", np.arange(400., 2401.)))[:bands]
        par = np.where(wavelengths <= 700., wavelengths * tr.PHOTONS, 0.)
        result["directPrad_leaf"] = float(direct[:, 2] @ par)
        for index, suffix in enumerate(("sunlit", "shaded")):
            result[f"diffusePrad_leaf_{suffix}"] = float(courtyard_diffuse[:, index + 2] @ par)
        if canopy is not None:
            nodes, _, _, leaf_area, _ = canopy
            result["diffusePrad_leaf"] = float((diffuse_node[:, ~nodes.soil].sum(axis=1) / leaf_area.sum()) @ par)
        else:
            result["diffusePrad_leaf"] = 0.
    _shares(result, geo)
    return result
