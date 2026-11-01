"""Four-temperature energy closure for slopes and statistical conical terrain.

Single slopes use horizontal-to-local beam projection and isotropic sky view.
Compound scenes use twelve cone azimuth populations plus the exposed plane,
physical surface areas, statistical mutual shadows and reciprocal lumped view
factors. The four temperatures remain shared by ALL facet populations. This is
an area-integrated approximation, not a spatially resolved terrain temperature
solver. Broadband inter-facet shortwave reflections use the supplied diffuse
spectral shape; view factors are conservative population factors, not a DEM
ray-tracing solution. LAI is per slope area for a single slope and per horizontal
area for compound terrain (facet LAI = horizontal LAI * cos(slope)).
"""
import math
from functools import lru_cache

import numpy as np

from .radiation import _shortwave_diffuse_fraction
from .solar import slope_direction
from .terrain_exchange import _hemisphere_sums, _mount_projection

SIGMA = 5.67e-8
_NAMES = ('soil_sunlit', 'soil_shaded', 'leaf_sunlit', 'leaf_shaded')


@lru_cache(maxsize=256)
def _geometry(kind, slope, aspect, height, radius, density, sza, saa):
    cosine = max(math.cos(math.radians(sza)), 0.)
    compound = kind == 2 and density > 0. and height > 0.
    if compound:
        occupied = density*math.pi*radius*radius
        if occupied >= 1.:
            raise ValueError('compound terrain requires density*pi*radius**2 < 1')
        beta = math.degrees(math.atan2(height, radius))
        cb = math.cos(math.radians(beta))
        area = np.array([1-occupied]+[occupied/cb/12]*12)
        slopes = (0.,)+(beta,)*12
        aspects = (0.,)+tuple(np.arange(12)*30.)
        tangent = math.tan(math.radians(min(max(sza, 0.), 89.99)))
        shadow = float(_mount_projection(height, radius, density, np.array([tangent]), outside=True)[0])
        visibility = np.array([math.exp(-shadow/(1-occupied))]+[
            -math.expm1(-shadow)/shadow if shadow > 1e-12 else 1.]*12)
        # A geometric seed from existing cone hemisphere transport. Symmetric
        # exchange areas then enforce reciprocity and exactly one horizontal
        # sky aperture, rather than independently clipping view-factor rows.
        mu, weights = np.polynomial.legendre.leggauss(12)
        mu, weights = (mu+1)/2, weights/2
        projection = _mount_projection(height, radius, density, np.sqrt(1-mu*mu)/mu, outside=True)
        ground_sky = float(np.sum(2*mu*weights*np.exp(-projection/(1-occupied))))
        upper, _ = _hemisphere_sums(np.array([[height, radius, density]]), 3)
        facet_sky = (1+cb)/2*float(np.mean(upper[0]))
        blocked = area*(1-np.array([ground_sky]+[facet_sky]*12))
        exchange = np.outer(blocked, blocked)
        np.fill_diagonal(exchange, 0.)
        target = float(area.sum()-1.)
        exchange *= target/max(float(exchange.sum()), 1e-300)
        if np.any(exchange.sum(axis=1) > area):
            # Positive conservative fallback for extreme cone aspect ratios.
            exchange = np.outer(area, area)
            np.fill_diagonal(exchange, 0.)
            exchange *= target/float(exchange.sum())
        factors = exchange/area[:, None]
        sky = 1-factors.sum(axis=1)
        if np.min(sky) < -1e-10:
            raise ValueError('terrain population geometry has no admissible sky view')
    else:
        slope = slope if kind == 1 else 0.
        aspect = aspect if kind == 1 else 0.
        area, slopes, aspects = np.array([1.]), (slope,), (aspect,)
        visibility = np.ones(1)
        sky = np.array([(1+math.cos(math.radians(slope)))/2])
        factors = np.zeros((1, 1))
    angles = [slope_direction(sza, saa, beta, azimuth)
              for beta, azimuth in zip(slopes, aspects)]
    local_z = np.array([angle[0] for angle in angles])
    local_a = np.array([angle[1] for angle in angles])
    visible = visibility*np.array([angle[2] and sza < 90. for angle in angles])
    direct = np.maximum(np.cos(np.radians(local_z)), 0.)*visible/max(cosine, 1e-8)
    if cosine <= 0.:
        direct[:] = 0.
    elif compound:
        # Every beam entering the horizontal aperture has exactly one first
        # intercept. Statistical shading determines its distribution, never
        # a duplicated projected-area source of extra sunlight.
        direct /= max(float(area @ direct), 1e-12)
    return dict(area=area, slopes=slopes, aspects=aspects, local_z=local_z,
                local_a=local_a, visibility=visible, direct=direct, sky=sky,
                factors=factors, compound=compound)


def layout(settings, global_sza, global_saa):
    """Cached immutable-by-convention geometry; arrays are never modified by callers."""
    terrain_type = 1 if settings.get('surface_model') == 'slope' else int(settings.get('terrain_type', 2))
    return _geometry(terrain_type,
                     float(settings.get('terrain_slope', settings.get('slope', 0.))),
                     float(settings.get('terrain_aspect', settings.get('aspect', 0.))),
                     max(float(settings.get('terrain_height', 20.)), 0.),
                     max(float(settings.get('terrain_radius', 10.)), 1e-6),
                     max(float(settings.get('terrain_density', 0.001)), 0.),
                     float(global_sza), float(global_saa))


def _patches(settings, geometry, plane_weights):
    patches = []
    for index, (beta, aspect) in enumerate(zip(geometry['slopes'], geometry['aspects'])):
        local = dict(settings, surface_model='plane')
        if geometry['compound']:
            cb = math.cos(math.radians(beta))
            local['lai'] = float(settings.get('lai', 0.))*cb
            if 'stand_density' in local:
                local['stand_density'] = float(local['stand_density'])*cb
        if str(local.get('vegetation_model', '')).lower() in ('row', 'raw'):
            local['row_azimuth'] = slope_direction(90., float(settings.get('row_azimuth', 0.)), beta, aspect)[1]
        z, az = min(float(geometry['local_z'][index]), 89.), float(geometry['local_a'][index])
        vis = float(geometry['visibility'][index])
        full = np.asarray(plane_weights(local, z, az, sunlit=vis > 0.), float)
        dark = np.array([0., full[0]+full[1], 0., full[2]+full[3]])
        patches.append((local, z, az, vis, full, dark))
    return patches


def flux_weights(settings, global_sza, global_saa, *, plane_weights):
    geometry = layout(settings, global_sza, global_saa)
    patches = _patches(settings, geometry, plane_weights)
    return tuple(sum(area*(vis*full+(1-vis)*dark)
                     for area, (_, _, _, vis, full, dark) in zip(geometry['area'], patches)))


def _scaled_sw(patch, meteo, direct, diffuse, callback):
    local, z, az, vis, _, _ = patch
    incoming = max(float(direct+diffuse), 0.)
    # Evaluate a common normalized linear spectral response and rescale it to
    # the actual facet input, including arbitrarily small positive irradiance.
    fraction = diffuse/incoming if incoming > 0. else 1.
    # This is a normalized linear response, not another real solar forcing.
    result = callback(dict(local, shortwave_diffuse_fraction=fraction,shortwave_enforce_solar_limit=0),
                      dict(meteo, shortwave=1000.), z, az, sunlit=vis > 0.)
    scale = incoming/1000.
    return {key: ([float(v)*scale for v in value] if isinstance(value, (list, tuple, np.ndarray))
                  else float(value) if key in ('shortwave_albedo', 'shortwave_diffuse_fraction')
                  else float(value)*scale)
            for key, value in result.items()}


def _sw_absorption(result, weights):
    return float(weights @ np.array([
        result.get('directVrad_soil', 0.)+result.get('diffuseVrad_soil_sunlit', 0.),
        result.get('diffuseVrad_soil_shaded', 0.),
        result.get('directVrad_leaf', 0.)+result.get('diffuseVrad_leaf_sunlit', 0.),
        result.get('diffuseVrad_leaf_shaded', 0.)]))


def _add_diffuse(base, unit, extra):
    result = dict(base)
    for key, value in unit.items():
        if key in ('shortwave_albedo', 'shortwave_diffuse_fraction'):
            continue
        if isinstance(value, list):
            result[key] = (np.asarray(base[key])+extra*np.asarray(value)).tolist()
        else:
            result[key] = base[key]+extra*value
    return result


def _densities(power, weights):
    density = np.divide(power, weights, out=np.zeros(4), where=weights > 1e-12)
    # A zero-area class contributes no energy. Give it the material's mean
    # diffuse input, matching the plane closure's continuous absent-class limit.
    for start in (0, 2):
        total_area = float(weights[start:start+2].sum())
        if total_area > 1e-12:
            mean = float(power[start:start+2].sum()/total_area)
            for index in (start, start+1):
                if weights[index] <= 1e-12:
                    density[index] = mean
    return density


def prepare_shortwave(settings, meteo, global_sza, global_saa, *, plane_shortwave, plane_weights):
    geometry = layout(settings, global_sza, global_saa)
    patches = _patches(settings, geometry, plane_weights)
    areas, sky, factors = (geometry[key] for key in ('area', 'sky', 'factors'))
    rin = max(float(meteo['shortwave']), 0.)
    fraction = _shortwave_diffuse_fraction(settings, rin, global_sza)
    horizontal_direct = rin*(1-fraction) if global_sza < 90. else 0.
    direct = horizontal_direct*geometry['direct']
    diffuse = rin*fraction*sky
    if not geometry['compound'] and sky[0] < 1.:
        ambient = dict(settings, surface_model='plane')
        ambient_z = min(float(global_sza), 89.)
        ambient_weights = np.asarray(plane_weights(ambient, ambient_z, global_saa, sunlit=global_sza < 90.), float)
        ambient_patch = (ambient, ambient_z, global_saa, float(global_sza < 90.), ambient_weights, ambient_weights)
        ambient_sw = _scaled_sw(ambient_patch, meteo, horizontal_direct, rin*fraction, plane_shortwave)
        diffuse += (1-sky)*max(rin-_sw_absorption(ambient_sw, ambient_weights), 0.)
    base = [_scaled_sw(patch, meteo, d, s, plane_shortwave) for patch, d, s in zip(patches, direct, diffuse)]
    if geometry['compound']:
        unit = [_scaled_sw(patch, meteo, 0., 1., plane_shortwave) for patch in patches]
        reflectivity = np.array([1-_sw_absorption(response, patch[4]) for response, patch in zip(unit, patches)])
        first_out = np.array([d+s-_sw_absorption(response, patch[4])
                              for d, s, response, patch in zip(direct, diffuse, base, patches)])
        radiosity = np.linalg.solve(np.eye(len(areas))-reflectivity[:, None]*factors, first_out)
        extra = factors @ radiosity
        base = [_add_diffuse(response, response_unit, added) for response, response_unit, added in zip(base, unit, extra)]
    weights = np.asarray(flux_weights(settings, global_sza, global_saa, plane_weights=plane_weights))
    direct_power, diffuse_power, par_direct, par_diffuse = np.zeros(4), np.zeros(4), 0., np.zeros(4)
    for area, response, (_, _, _, vis, full, _) in zip(areas, base, patches):
        for a, component in ((0, 'soil'), (2, 'leaf')):
            direct_power[a] += area*full[a]*response.get('directVrad_'+component, 0.)
            ds = response.get('diffuseVrad_'+component+'_sunlit', 0.)
            dh = response.get('diffuseVrad_'+component+'_shaded', 0.)
            diffuse_power[a] += area*vis*full[a]*ds
            diffuse_power[a+1] += area*(full[a+1]*dh+(1-vis)*full[a]*ds)
        par_direct += area*full[2]*response.get('directPrad_leaf', 0.)
        ps = response.get('diffusePrad_leaf_sunlit', response.get('diffusePrad_leaf', 0.))
        ph = response.get('diffusePrad_leaf_shaded', response.get('diffusePrad_leaf', 0.))
        par_diffuse[2] += area*vis*full[2]*ps
        par_diffuse[3] += area*(full[3]*ph+(1-vis)*full[2]*ps)
    direct_density = np.divide(direct_power, weights, out=np.zeros(4), where=weights > 1e-12)
    diffuse_density, par_density = _densities(diffuse_power, weights), _densities(par_diffuse, weights)
    result = {}
    for a, component in ((0, 'soil'), (2, 'leaf')):
        result['directVrad_'+component] = float(direct_density[a])
        result['diffuseVrad_'+component+'_sunlit'] = float(diffuse_density[a])
        result['diffuseVrad_'+component+'_shaded'] = float(diffuse_density[a+1])
        result['diffuseVrad_'+component] = float(diffuse_power[a:a+2].sum()/max(weights[a:a+2].sum(), 1e-12))
    result.update(directPrad_leaf=float(par_direct/max(weights[2], 1e-12)),
                  diffusePrad_leaf_sunlit=float(par_density[2]), diffusePrad_leaf_shaded=float(par_density[3]),
                  diffusePrad_leaf=float(par_diffuse.sum()/max(weights[2:].sum(), 1e-12)))
    # Spectrum fields retain their original per-soil-area convention. Incident
    # spectra are scene powers per ground area; they exclude no illuminated facet.
    for key in ('directVrad_incident_spectrum', 'diffuseVrad_incident_spectrum',
                'diffuseVrad_soil_spectrum', 'diffuseVrad_soil_incident_spectrum'):
        spectra = [np.asarray(response.get(key, []), float) for response in base]
        result[key] = sum(area*spec for area, spec in zip(areas, spectra)).tolist()
        if 'soil' in key:
            result[key] = (np.asarray(result[key])/max(weights[:2].sum(), 1e-12)).tolist()
    direct_incident = float(areas @ direct)
    external_diffuse = float(areas @ (rin*fraction*sky))
    for key in ('directVrad_incident', 'diffuseVrad_incident'):
        result[key] = sum(float(area)*response.get(key, 0.) for area, response in zip(areas, base))
    absorbed = float(direct_power.sum()+diffuse_power.sum())
    external = direct_incident+(external_diffuse if geometry['compound'] else float(areas @ diffuse))
    result.update(shortwave_global=rin, shortwave_direct_incident=direct_incident,
                  shortwave_diffuse_incident=(external_diffuse if geometry['compound'] else float(areas @ diffuse)),
                  shortwave_diffuse_fraction=fraction,
                  shortwave_albedo=1-absorbed/max(external, 1e-12) if external > 0. else 0.,
                  terrain_direct_intercepted_ground=direct_incident,
                  terrain_diffuse_sky_ground=external_diffuse,
                  terrain_soil_area_ground=float(weights[:2].sum()),
                  terrain_flux_weights=tuple(float(value) for value in weights))
    return result


def _lw_response(patch, meteo, state, irradiance, callback):
    local, z, az, vis, full, dark = patch
    result = np.zeros(4)
    if vis > 0.:
        lit = callback(local, dict(meteo, longwave=irradiance), state, z, az, sunlit=True)
        result += vis*full*np.array([lit['LW_'+name] for name in _NAMES])
    if vis < 1.:
        shade = callback(local, dict(meteo, longwave=irradiance), state, z, az, sunlit=False)
        result += (1-vis)*dark*np.array([shade['LW_'+name] for name in _NAMES])
    return result


def _lw_source(patch, meteo, state, callback):
    p0 = _lw_response(patch, meteo, state, 0., callback)
    p1 = _lw_response(patch, meteo, state, 1., callback)-p0
    local, _, _, vis, full, dark = patch
    weights = vis*full+(1-vis)*dark
    es, el = float(local['soil_emissivity']), float(local['leaf_emissivity'])
    emitted = weights @ (SIGMA*np.array([es*state['Tss']**4, es*state['Tsh']**4,
                                         2*el*state['Tls']**4, 2*el*state['Tlh']**4]))
    return p0, p1, float(emitted-p0.sum()), float(1-p1.sum())


def prepare_longwave(settings, meteo, state, global_sza, global_saa, *, plane_longwave, plane_weights):
    geometry = layout(settings, global_sza, global_saa)
    patches = _patches(settings, geometry, plane_weights)
    areas, sky, factors = (geometry[key] for key in ('area', 'sky', 'factors'))
    sky_lw = max(float(meteo['longwave']), 0.)
    responses = [_lw_source(patch, meteo, state, plane_longwave) for patch in patches]
    source = np.array([response[2] for response in responses])
    reflection = np.array([response[3] for response in responses])
    if geometry['compound']:
        radiosity = np.linalg.solve(np.eye(len(areas))-reflection[:, None]*factors,
                                    source+reflection*sky*sky_lw)
        incident = sky*sky_lw+factors @ radiosity
    else:
        ambient_geometry = dict(geometry, slopes=(0.,), aspects=(0.,), local_z=np.array([min(global_sza, 89.)]),
                                local_a=np.array([global_saa]), visibility=np.array([float(global_sza < 90.)]), compound=False)
        ambient = _patches(settings, ambient_geometry, plane_weights)[0]
        _, _, ambient_source, ambient_reflection = _lw_source(ambient, meteo, state, plane_longwave)
        incident = sky*sky_lw+(1-sky)*(ambient_source+ambient_reflection*sky_lw)
    power = sum(area*(p0+p1*incoming) for area, (p0, p1, _, _), incoming in zip(areas, responses, incident))
    weights = np.asarray(flux_weights(settings, global_sza, global_saa, plane_weights=plane_weights))
    density = _densities(power, weights)
    result = {'LW_'+name: float(value) for name, value in zip(_NAMES, density)}
    result.update(longwave_down=sky_lw, tir_scatter_soil=float(density[:2].sum()/2),
                  tir_scatter_leaf=float(density[2:].sum()/2))
    return result
