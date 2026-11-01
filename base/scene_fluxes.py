"""Scene budgets per declared reference area.

Plane/compound terrain/urban use horizontal ground area; a single slope uses
its local slope-plane area consistently with its LAI and incident irradiance.
Monin-Obukhov stability needs fluxes per unit horizontal area, so a single
slope's budgets are divided by cos(slope) for that purpose only.
"""
import math

from .tirteb_radiation import flux_weights


def horizontal_area_factor(settings):
    """Reference area per unit horizontal area: 1/cos(slope) for a single slope."""
    terrain_type = 1 if settings.get('surface_model') == 'slope' else int(settings.get('terrain_type', 2))
    if settings.get('surface_model') in ('terrain', 'slope') and terrain_type == 1:
        slope = float(settings.get('terrain_slope', settings.get('slope', 0.0)))
        return 1.0 / max(math.cos(math.radians(slope)), 1e-6)
    return 1.0


def urban_weights(settings, state):
    density = float(settings.get('urban_density', 0.0015))
    length = float(settings.get('urban_building_length', 10.0))
    width = float(settings.get('urban_building_width', 10.0))
    height = float(settings.get('urban_building_height', 20.0))
    roof = density * length * width
    wall = density * 2.0 * height * (length + width)
    ground = 1.0 - roof
    if ground <= 0.0:
        raise ValueError('Building footprint must be less than one.')
    lit = float(state.get('sza', 90.0)) < 90.0 and state.get('sunlit_surface', True)
    # Incoming radiation stores the area-mean building shadow shares.
    roof_sun = float(state.get('roof_sunlit_fraction', float(lit)))
    wall_sun = float(state.get('wall_sunlit_fraction', 0.5 if lit else 0.0))
    building_sun = float(state.get('street_building_sunlit_fraction', float(lit)))
    ps = float(state.get('soil_sunlit_fraction', 1.0 if lit else 0.0)) * building_sun
    lai = float(settings.get('lai', 0.0))
    leaf_sun = float(state.get('lai_sunlit', 0.0)) / lai if lai > 0 else 0.0
    leaf_shade = float(state.get('lai_shaded', 0.0)) / lai if lai > 0 else 0.0
    return {'roofsunlit': roof * roof_sun, 'roofshaded': roof * (1.0-roof_sun),
            'wallsunlit': wall * wall_sun, 'wallshaded': wall * (1.0-wall_sun),
            'streetsunlit': ground * ps, 'streetshaded': ground * (1.0-ps),
            'leafsunlit': ground * leaf_sun * building_sun,
            'leafshaded': ground * (leaf_shade + leaf_sun * (1.0-building_sun))}


def component_weights(settings, state, urban=False):
    if urban:
        return urban_weights(settings, state)
    soil_sun, soil_shade, leaf_sun, leaf_shade = flux_weights(
        settings, state.get('sza_global',state['sza']),
        state.get('saa_global',state['saa']),
        float(state.get('sza_global',state['sza'])) < 90.0)
    lai = float(settings.get('lai', 0.0))
    return dict(zip(('soilsunlit','soilshaded','leafsunlit','leafshaded'),
                    (soil_sun, soil_shade, leaf_sun/lai if lai else 0.0,
                     leaf_shade/lai if lai else 0.0)))


def scene_totals(settings, state, urban=False):
    weights = component_weights(settings, state, urban)
    totals = {quantity: sum(weight * float(state.get(prefix + name, 0.0))
                           for name, weight in weights.items())
              for quantity, prefix in (('Rn','N_'),('LE','LE_'),('H','H_'))}
    totals['G'] = sum(weight * float(state.get(
        ('G_sunlit' if name == 'soilsunlit' else 'G_shaded') if name.startswith('soil')
        else 'G_' + name, 0.0)) for name, weight in weights.items())
    return totals
