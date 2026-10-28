"""TIRTEB radiation budget and observations with the TIRT standard canopy closure.

Every canopy quantity comes from canopy_closure (standard four components) and node_radiosity
(hom3 / row3 / crown2 schedules): view fractions, the beam reaching each node, and all orders of
shortwave and longwave exchange between sky, canopy nodes and soil nodes.

Units of the component budget: soil per unit soil area; leaves per unit leaf area. The energy balance
treats the canopy as one big leaf per unit ground (bulk resistances), so it multiplies the leaf terms by
LAI; scene fluxes weight soil by Ps / 1 - Ps (sunlit soil seen from the sun) and leaves by
LAI_sunlit / LAI and LAI_shaded / LAI, with LAI_sunlit = intercepted beam / (G / cos sza).
"""
import functools
import math

import numpy as np

from . import canopy_closure, node_radiosity
from .radiation import (_shortwave_diffuse_fraction, _requested_shortwave_diffuse_fraction,
                        _solar_direct_normal_limit, clearness_index, optical_properties)
from .solar import slope_direction

G = 0.5
CLOSURE = "standard"
SIGMA = 5.67e-8
BANDS = 2001
# W m-2 per band (1 nm) * wavelength (nm) -> umol photons m-2 s-1.
PHOTONS = 1e-3 / (6.02214e23 * 6.6262e-34 * 299792458.0)
_KINDS = {"homogeneous": "hom", "raw": "row", "forest": "crown", "bare_soil": "bare", "soil": "bare", "barren": "bare"}


def _slope(settings):
    """(slope, aspect) of a single-slope scene, else None."""
    if settings.get("surface_model") == "terrain" and int(settings.get("terrain_type", 1)) == 1:
        return float(settings.get("terrain_slope", 0.0)), float(settings.get("terrain_aspect", 0.0))
    return None


def canopy(settings):
    """(kind, lai, hspot, shape, row azimuth in the surface frame, crown base) of the TIRTEB canopy."""
    kind = str(settings.get("vegetation_model", "hom")).strip().lower()
    kind = _KINDS.get(kind, kind)
    lai, hspot = max(float(settings.get("lai", 0.0)), 0.0), float(settings.get("hspot", 0.1))
    if kind == "crown":
        shape = (float(settings.get("stand_density", 0.02)), float(settings.get("height", 3.0)),
                 float(settings.get("width", 3.0)))
        base = settings.get("crown_base_height")
        if shape[0] > 0.0 and lai > 0.0:
            return kind, lai, hspot, shape, 0.0, shape[1] if base in (None, "") else float(base)
    elif kind == "row" and lai > 0.0:
        shape = tuple(float(settings.get(k, d)) for k, d in (("row_width", 0.5), ("row_blank", 0.5), ("row_height", 1.0)))
        azimuth = float(settings.get("row_azimuth", 0.0))
        slope = _slope(settings)
        if slope is not None:
            azimuth = float(slope_direction(90.0, azimuth, *slope)[1])
        return kind, lai, hspot, shape, azimuth, None
    elif kind == "hom" and lai > 0.0:
        return kind, lai, hspot, canopy_closure.HOMOGENEOUS, 0.0, None
    return "bare", 0.0, hspot, canopy_closure.HOMOGENEOUS, 0.0, None


def _shade(parts):
    return np.array([0.0, parts[0] + parts[1], 0.0, parts[2] + parts[3]])


def _plane_fractions(settings, vza, vaa, sza, saa, sunlit=True):
    kind, lai, hspot, shape, azimuth, base = canopy(settings)
    if kind == "bare":
        parts = np.array([1.0, 0.0, 0.0, 0.0])
    else:
        parts = canopy_closure.fractions(kind, lai, hspot, shape, vza, vaa, min(float(sza), 89.9), saa, azimuth,
                                         base, CLOSURE)
    return parts if sunlit else _shade(parts)


def view_fractions(settings, vza, vaa, sza, saa, sunlit=True):
    """Standard four components [soil sun, soil shade, leaf sun, leaf shade] along a view (surface frame)."""
    if settings.get("surface_model") == "terrain" and int(settings.get("terrain_type", 1)) == 2:
        from .terrain_rt import compound_terrain_proportions
        return compound_terrain_proportions(
            settings, sza, saa, vza, vaa,
            lambda s, sa, v, va, lit: _plane_fractions(settings, v, va, s, sa, lit), sunlit)
    return _plane_fractions(settings, vza, vaa, sza, saa, sunlit)


@functools.lru_cache(maxsize=4096)
def _layout(kind, lai, shape, sza, saa, azimuth, sunlit):
    """Nodes, unit beam reaching each node, soil and leaf area per node, sunlit area per node."""
    nodes, beam = node_radiosity.canopy_nodes(canopy_closure.NODES[kind], lai, shape, min(sza, 89.9), saa, azimuth)
    beam = np.asarray(beam, float) if sunlit and sza < 90.0 else np.zeros(len(nodes.capacity))
    extinction = G / math.cos(math.radians(min(sza, 85.0)))
    soil_area = np.where(nodes.soil, nodes.capacity, 0.0)
    leaf_area = np.where(nodes.soil, 0.0, nodes.capacity / (4.0 * G))
    sunlit_area = np.where(nodes.soil, np.minimum(beam, soil_area), np.minimum(beam / extinction, leaf_area))
    return nodes, beam, soil_area, leaf_area, sunlit_area


def _scene(settings, sza, saa, sunlit):
    kind, lai, _, shape, azimuth, _ = canopy(settings)
    if kind == "bare":
        return None
    return _layout(kind, lai, tuple(shape), float(sza), float(saa), azimuth, bool(sunlit and float(sza) < 90.0))


def _plane_flux_weights(settings, sza, saa, sunlit=True):
    """Per unit ground: sunlit and shaded soil fractions and sunlit and shaded leaf area indices."""
    layout = _scene(settings, sza, saa, sunlit)
    lit = bool(sunlit) and float(sza) < 90.0
    if layout is None:
        return (1.0, 0.0, 0.0, 0.0) if lit else (0.0, 1.0, 0.0, 0.0)
    nodes, _, _, leaf_area, sunlit_area = layout
    soil = float(sunlit_area[nodes.soil].sum())
    leaf = float(sunlit_area[~nodes.soil].sum())
    lai = float(leaf_area.sum())
    return soil, 1.0 - soil, leaf, lai - leaf


def _means(per_node, area, sunlit_area, mask):
    """Absorption per unit area of the sunlit and shaded parts of the masked nodes (bands x 2)."""
    density = np.divide(per_node, area, out=np.zeros_like(per_node), where=area > 0.0)
    result = []
    for weights in (np.where(mask, sunlit_area, 0.0), np.where(mask, area - sunlit_area, 0.0)):
        if weights.sum() <= 1.0e-12:
            weights = np.where(mask, area, 0.0)
        result.append(density @ (weights / weights.sum()))
    return np.stack(result, axis=-1)


def _spectrum(settings, key, fallback):
    values = settings.get(key, [])
    if len(values) >= BANDS:
        return np.asarray(values[:BANDS], float)
    return np.full(BANDS, float(fallback))


def _empty_shortwave():
    result = {key: 0.0 for key in (
        "directVrad_incident", "diffuseVrad_incident", "directVrad_leaf", "directVrad_soil",
        "diffuseVrad_leaf", "diffuseVrad_soil", "diffuseVrad_leaf_sunlit", "diffuseVrad_leaf_shaded",
        "diffuseVrad_soil_sunlit", "diffuseVrad_soil_shaded", "directPrad_leaf", "diffusePrad_leaf",
        "shortwave_global", "shortwave_direct_incident", "shortwave_diffuse_incident", "shortwave_diffuse_fraction",
        "shortwave_albedo", "diffusePrad_leaf_sunlit", "diffusePrad_leaf_shaded")}
    for key in ("directVrad_incident_spectrum", "diffuseVrad_incident_spectrum", "diffuseVrad_soil_spectrum",
                "diffuseVrad_soil_incident_spectrum"):
        result[key] = []
    return result


SHORTWAVE_BIN_NM = 50


def _spectral_groups(bin_nm=SHORTWAVE_BIN_NM):
    """Contiguous 1 nm band groups used by the energy-balance shortwave solver."""
    width = max(1, int(bin_nm))
    return [(start, min(start + width, BANDS)) for start in range(0, BANDS, width)]


def _group_sum(values, groups):
    return np.asarray([float(np.sum(values[start:end])) for start, end in groups], float)


def _group_mean(values, groups):
    return np.asarray([float(np.mean(values[start:end])) for start, end in groups], float)


def _expand_groups(values, groups):
    """Expand group totals back to the 1 nm grid without changing each group total."""
    expanded = np.zeros(BANDS, float)
    for value, (start, end) in zip(values, groups):
        expanded[start:end] = float(value) / (end - start)
    return expanded


def _plane_prepare_shortwave(settings, meteo, sza, saa, sunlit=True):
    """Absorbed shortwave: beam per unit sunlit area (soil, leaf) and diffuse per unit area of every component."""
    result = _empty_shortwave()
    rin = float(meteo["shortwave"])
    if rin <= 0.0:
        return result
    fraction = _shortwave_diffuse_fraction(settings, rin, sza)
    result.update(shortwave_global=rin, shortwave_direct_incident=rin * (1.0 - fraction),
                  shortwave_diffuse_incident=rin * fraction, shortwave_diffuse_fraction=fraction)
    rho, tau, soil, _, _ = optical_properties(settings)
    fesun, fesky = (np.asarray(settings.get(key, [])[:BANDS], float) for key in ("fesun", "fesky"))
    wavelengths = np.asarray(settings.get("wavelengths", [])[:BANDS], float)
    if fesun.size != BANDS or fesky.size != BANDS:
        # Uniform 1 nm bands; the 0.001 below converts the per-um weights of the solar files.
        fesun = fesky = np.full(BANDS, 1000.0 / BANDS)
        wavelengths = np.arange(400.0, 2401.0)
    lit = bool(sunlit)
    groups = _spectral_groups()
    esun_full = rin * (1.0 - fraction) * fesun * 0.001
    esky_full = rin * fraction * fesky * 0.001
    esun = _group_sum(esun_full, groups)
    esky = _group_sum(esky_full, groups)
    wavelengths = _group_mean(wavelengths, groups)
    beam = esun if lit else np.zeros(len(groups))
    leaf_rho_full = np.clip(_spectrum(settings, "leaf_reflectance_spectrum", rho), 0.0, 1.0)
    leaf_tau_full = np.minimum(
        np.clip(_spectrum(settings, "leaf_transmittance_spectrum", tau), 0.0, 1.0), 1.0 - leaf_rho_full
    )
    leaf_rho = _group_mean(leaf_rho_full, groups)
    leaf_tau = _group_mean(leaf_tau_full, groups)
    background = ("urban_reflectance_street_spectrum"
                  if str(settings.get("surface_model", "")).strip().lower() == "urban" else "soil_reflectance_spectrum")
    soil_rho = _group_mean(np.clip(_spectrum(settings, background, soil), 0.0, 1.0), groups)
    result.update(directVrad_incident=float(beam.sum()),
                  directVrad_incident_spectrum=_expand_groups(beam, groups).tolist(),
                  diffuseVrad_incident=float(esky.sum()),
                  diffuseVrad_incident_spectrum=_expand_groups(esky, groups).tolist())
    par = np.where((wavelengths >= 400.0) & (wavelengths <= 700.0), wavelengths * PHOTONS, 0.0)
    layout = _scene(settings, sza, saa, lit)
    if layout is None:
        diffuse = esky * (1.0 - soil_rho)
        result.update(directVrad_soil=float((beam * (1.0 - soil_rho)).sum()), diffuseVrad_soil=float(diffuse.sum()),
                      diffuseVrad_soil_sunlit=float(diffuse.sum()), diffuseVrad_soil_shaded=float(diffuse.sum()),
                      diffuseVrad_soil_spectrum=_expand_groups(diffuse, groups).tolist(),
                      diffuseVrad_soil_incident_spectrum=_expand_groups(esky, groups).tolist(),
                      shortwave_albedo=float(((beam + esky) * soil_rho).sum() / max((beam + esky).sum(), 1e-12)))
        return result
    nodes, unit_beam, soil_area, leaf_area, sunlit_area = layout
    beam_first, beam_scattered = node_radiosity.absorbed(nodes, soil_rho, leaf_rho, leaf_tau,
                                                         node_radiosity.from_above(nodes, unit_beam))
    sky_first, sky_scattered = node_radiosity.absorbed(nodes, soil_rho, leaf_rho, leaf_tau,
                                                       node_radiosity.sky_arrivals(nodes))
    direct = beam[:, None] * beam_first
    diffuse = beam[:, None] * beam_scattered + esky[:, None] * (sky_first + sky_scattered)
    soil_sunlit, leaf_sunlit = (float(sunlit_area[mask].sum()) for mask in (nodes.soil, ~nodes.soil))
    direct_soil = direct[:, nodes.soil].sum(axis=1) / soil_sunlit if soil_sunlit > 1e-12 else np.zeros(len(groups))
    direct_leaf = direct[:, ~nodes.soil].sum(axis=1) / leaf_sunlit if leaf_sunlit > 1e-12 else np.zeros(len(groups))
    soil_parts = _means(diffuse, soil_area, sunlit_area, nodes.soil)
    leaf_parts = _means(diffuse, leaf_area, sunlit_area, ~nodes.soil)
    soil_mean = diffuse[:, nodes.soil].sum(axis=1) / soil_area.sum()
    leaf_mean = diffuse[:, ~nodes.soil].sum(axis=1) / leaf_area.sum()
    absorbed_total = float((direct + diffuse).sum())
    soil_incident = soil_mean / np.maximum(1.0 - soil_rho, 1e-12)
    result.update(
        directVrad_soil=float(direct_soil.sum()), directVrad_leaf=float(direct_leaf.sum()),
        diffuseVrad_soil=float(soil_mean.sum()), diffuseVrad_leaf=float(leaf_mean.sum()),
        diffuseVrad_soil_sunlit=float(soil_parts[:, 0].sum()), diffuseVrad_soil_shaded=float(soil_parts[:, 1].sum()),
        diffuseVrad_leaf_sunlit=float(leaf_parts[:, 0].sum()), diffuseVrad_leaf_shaded=float(leaf_parts[:, 1].sum()),
        diffuseVrad_soil_spectrum=_expand_groups(soil_mean, groups).tolist(),
        diffuseVrad_soil_incident_spectrum=_expand_groups(soil_incident, groups).tolist(),
        directPrad_leaf=float((direct_leaf * par).sum()), diffusePrad_leaf=float((leaf_mean * par).sum()),
        diffusePrad_leaf_sunlit=float((leaf_parts[:, 0] * par).sum()),
        diffusePrad_leaf_shaded=float((leaf_parts[:, 1] * par).sum()),
        shortwave_albedo=1.0 - absorbed_total / max(float((beam + esky).sum()), 1e-12),
    )
    return result


def longwave(settings, sky, temperatures, sza, saa, sunlit=True, soil_emissivity=None):
    """Absorbed longwave per unit area [soil sun, soil shade, leaf sun, leaf shade] from an isotropic downward
    irradiance `sky` and the emission of all nodes at temperatures (Tss, Tsh, Tls, Tlh), all orders."""
    es = max(float(settings["soil_emissivity"] if soil_emissivity is None else soil_emissivity), 0.0)
    el = max(float(settings["leaf_emissivity"]), 0.0)
    layout = _scene(settings, sza, saa, sunlit)
    if layout is None:
        return np.array([es * sky, es * sky, 0.0, 0.0])
    nodes, _, soil_area, leaf_area, sunlit_area = layout
    area = soil_area + leaf_area
    share = np.divide(sunlit_area, area, out=np.zeros_like(area), where=area > 0.0)
    tss, tsh, tls, tlh = (float(t) ** 4 for t in temperatures)
    exitance = np.where(nodes.soil, es * (share * tss + (1.0 - share) * tsh),
                        el * (share * tls + (1.0 - share) * tlh)) * SIGMA / math.pi
    first, scattered = node_radiosity.absorbed(nodes, 1.0 - es, 1.0 - el, 0.0,
                                               float(sky) * node_radiosity.sky_arrivals(nodes),
                                               nodes.emission(exitance[:, None])[:, 0])
    total = first + scattered
    soil = _means(total, soil_area, sunlit_area, nodes.soil)[0]
    leaf = _means(total, leaf_area, sunlit_area, ~nodes.soil)[0]
    return np.array([soil[0], soil[1], leaf[0], leaf[1]])


def _plane_prepare_longwave(settings, meteo, state, sza, saa, sunlit=True):
    sky = float(meteo["longwave"])
    values = longwave(settings, sky, (state["Tss"], state["Tsh"], state["Tls"], state["Tlh"]), sza, saa, sunlit)
    return {"LW_soil_sunlit": values[0], "LW_soil_shaded": values[1], "LW_leaf_sunlit": values[2],
            "LW_leaf_shaded": values[3], "longwave_down": sky,
            "tir_scatter_soil": 0.5 * (values[0] + values[1]), "tir_scatter_leaf": 0.5 * (values[2] + values[3])}


def flux_weights(settings, sza, saa, sunlit=True):
    """Physical component areas; terrain angles are in the global frame."""
    if settings.get('surface_model') in ('terrain','slope'):
        from .terrain_energy import flux_weights as terrain_weights
        return terrain_weights(settings, sza, saa, plane_weights=_plane_flux_weights)
    return _plane_flux_weights(settings, sza, saa, sunlit)


def prepare_shortwave(settings, meteo, sza, saa, sunlit=True):
    if settings.get('surface_model') in ('terrain','slope'):
        from .terrain_energy import prepare_shortwave as terrain_shortwave
        result = terrain_shortwave(settings, meteo, sza, saa,
                                 plane_shortwave=_plane_prepare_shortwave,
                                 plane_weights=_plane_flux_weights)
    else:
        result = _plane_prepare_shortwave(settings, meteo, sza, saa, sunlit)
    requested = _requested_shortwave_diffuse_fraction(settings,float(meteo['shortwave']),sza)
    actual = float(result.get('shortwave_diffuse_fraction',requested))
    kt = clearness_index(settings, float(meteo['shortwave']), sza)
    result.update(shortwave_diffuse_fraction_requested=requested,
                  shortwave_partition_limited=actual>requested+1e-12,
                  solar_direct_normal_limit_W_m2=_solar_direct_normal_limit(settings, sza),
                  clearness_index=kt, clearness_index_exceeds_one=bool(kt > 1.0),
                  shortwave_direct_normal_W_m2=(float(meteo['shortwave'])*(1-actual)/max(math.cos(math.radians(sza)),1e-12)
                                                if sza<90 and float(meteo['shortwave'])>0 else 0.0))
    return result


def prepare_longwave(settings, meteo, state, sza, saa, sunlit=True):
    if settings.get('surface_model') in ('terrain','slope'):
        from .terrain_energy import prepare_longwave as terrain_longwave
        return terrain_longwave(settings, meteo, state, sza, saa,
                                plane_longwave=_plane_prepare_longwave,
                                plane_weights=_plane_flux_weights)
    return _plane_prepare_longwave(settings, meteo, state, sza, saa, sunlit)


def merge_radiation(shortwave, longwave_terms):
    """Total absorbed radiation per unit area of each component."""
    result = dict(shortwave)
    result.update(longwave_terms)
    result["shortwave_direct"] = shortwave.get("directVrad_leaf", 0.0) + shortwave.get("directVrad_soil", 0.0)
    result["shortwave_diffuse"] = shortwave.get("diffuseVrad_leaf", 0.0) + shortwave.get("diffuseVrad_soil", 0.0)
    for component in ("soil", "leaf"):
        result[f"{component}_sunlit"] = (shortwave.get(f"directVrad_{component}", 0.0)
                                         + shortwave.get(f"diffuseVrad_{component}_sunlit", 0.0)
                                         + longwave_terms.get(f"LW_{component}_sunlit", 0.0))
        result[f"{component}_shaded"] = (shortwave.get(f"diffuseVrad_{component}_shaded", 0.0)
                                         + longwave_terms.get(f"LW_{component}_shaded", 0.0))
    return result


def diagnostics(settings, state, vza, vaa):
    """Flux weights, view/sun gaps, view four components and node-radiosity thermal scattering."""
    sza, saa, sunlit = float(state["sza"]), float(state["saa"]), bool(state.get("sunlit_surface", True))
    weights = flux_weights(settings, float(state.get('sza_global',sza)),
                           float(state.get('saa_global',saa)), sunlit)
    view = view_fractions(settings, vza, vaa, sza, saa, sunlit and sza < 90.0)
    kind, lai, _, shape, azimuth, _ = canopy(settings)
    scatter = np.zeros(4)
    if kind != "bare":
        scatter = node_radiosity.thermal_scatter(canopy_closure.NODES[kind], lai, shape,
                                                 float(settings["leaf_emissivity"]), float(settings["soil_emissivity"]),
                                                 vza, vaa, min(sza, 89.9), saa, azimuth)
        if not (sunlit and sza < 90.0):
            scatter = _shade(scatter)
    return {
        "soil_sunlit_fraction": weights[0], "soil_shaded_fraction": weights[1],
        "lai_sunlit": weights[2], "lai_shaded": weights[3],
        "gap_view": float(view[0] + view[1]), "gap_sun": weights[0],
        "tirt_direct_soil_sun": view[0], "tirt_direct_soil_shade": view[1],
        "tirt_direct_leaf_sun": view[2], "tirt_direct_leaf_shade": view[3],
        "tirt_scatter_soil_sun": scatter[0], "tirt_scatter_soil_shade": scatter[1],
        "tirt_scatter_leaf_sun": scatter[2], "tirt_scatter_leaf_shade": scatter[3],
    }


def reflectance(settings, sza, saa, vza, vaa, rho, tau, soil, diffuse=0.0, sunlit=True, soil_direct=None):
    """Plane or single-slope BRF / HDRF (surface-frame angles) of the standard closure."""
    kind, lai, hspot, shape, azimuth, base = canopy(settings)
    if kind == "bare":
        kind, lai = "hom", 0.0
    return canopy_closure.reflectance(kind, lai, hspot, shape, float(rho), float(tau), float(soil), float(vza),
                                      float(vaa), min(float(sza), 89.9), float(saa), azimuth, base, float(diffuse),
                                      bool(sunlit), None if soil_direct is None else float(soil_direct), CLOSURE)


def slope_reflectance(settings, global_sza, global_saa, sza, saa, vza, vaa, rho, tau, soil, diffuse=0.0,
                      sunlit=True, soil_direct=None):
    """Single-slope BRF referred to the horizontal global irradiance, as base.optical: the local direct
    and sky reflectances weighted by the energy-balance projection (cos i / cos sza, (1 + cos slope) / 2)."""
    from .terrain_energy import layout
    geometry = layout(settings, global_sza, global_saa)
    direct = (1.0 - float(diffuse)) * float(geometry["direct"][0]) if sunlit else 0.0
    sky = float(diffuse) * float(geometry["sky"][0])
    value = 0.0
    if direct > 0.0:
        value += direct * reflectance(settings, sza, saa, vza, vaa, rho, tau, soil, 0.0, True, soil_direct)
    if sky > 0.0:
        value += sky * reflectance(settings, sza, saa, vza, vaa, rho, tau, soil, 1.0, True, soil_direct)
    return value
