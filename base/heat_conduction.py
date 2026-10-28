"""Heat diffusion with explicit units and elapsed-time thermal memory.

Soil ``cs`` is specific heat [J kg-1 K-1]; urban ``cs`` is already volumetric
[J m-3 K-1]. Positive conductive flux enters the material from its surface.
"""
import math

HISTORY_SECONDS = 48.0 * 3600.0
# Soil/street column for the depth method: 2 m deep with a zero-flux bottom,
# about 11 diurnal damping depths for typical soil diffusivity (~1e-6 m2 s-1).
COLUMN_DEPTH = 2.0
COLUMN_TOP_LAYER = 0.01
COLUMN_LAYER_GROWTH = 1.3


def thermal_properties(settings, urban=False):
    """Return conductivity [W m-1 K-1], volumetric capacity [J m-3 K-1]."""
    prefix = "urban_" if urban else "soil_"
    conductivity = max(float(settings.get(prefix + "lambda", 1.5 if urban else 1.55)), 0.0)
    capacity = float(settings.get(prefix + "cs", 1000000.0 if urban else 1180.0))
    if not urban:
        capacity *= float(settings.get("soil_rhos", 1100.0))
    return conductivity, max(capacity, 1.0e-6)


def time_step(settings, state):
    return max(float(state.get("dt_seconds", settings.get("dtime", 1800.0))), 1.0)


def implicit_profile(old, surface, bottom, conductivity, capacity, dz, dt):
    """Backward Euler solve with current Dirichlet surface and bottom values.

    For interior node i: (1+2r) T_i - r T_(i-1) - r T_(i+1) = T_i_old,
    r = conductivity * dt / (capacity * dz**2). No Fourier-number clipping.
    The boundary nodes carry no storage; interior node storage has width dz.
    """
    old = [float(value) for value in old]
    count = len(old) - 2
    if count <= 0:
        return [float(surface), float(bottom)]
    ratio = conductivity * dt / (capacity * dz * dz)
    diagonal = [1.0 + 2.0 * ratio] * count
    rhs = old[1:-1]
    rhs[0] += ratio * float(surface)
    rhs[-1] += ratio * float(bottom)
    for index in range(1, count):
        multiplier = -ratio / diagonal[index - 1]
        diagonal[index] += multiplier * ratio
        rhs[index] -= multiplier * rhs[index - 1]
    result = [0.0] * count
    result[-1] = rhs[-1] / diagonal[-1]
    for index in range(count - 2, -1, -1):
        result[index] = (rhs[index] + ratio * result[index + 1]) / diagonal[index]
    return [float(surface)] + result + [float(bottom)]


def column_layers():
    """Layer thicknesses [m] of the soil/street column, growing geometrically with depth."""
    layers = []
    total = 0.0
    thickness = COLUMN_TOP_LAYER
    while total + thickness < COLUMN_DEPTH - 1.0e-9:
        layers.append(thickness)
        total += thickness
        thickness *= COLUMN_LAYER_GROWTH
    remainder = COLUMN_DEPTH - total
    if layers and remainder < layers[-1]:
        layers[-1] += remainder
    else:
        layers.append(remainder)
    return layers


def column_profile(old, surface, conductivity, capacity, layers, dt):
    """Backward Euler finite volumes: Dirichlet surface, zero-flux bottom.

    Node i stores capacity*(h_(i-1)+h_i)/2; the bottom node stores half its layer.
    """
    count = len(layers)
    conductance = [conductivity / h for h in layers]
    storage = [capacity * 0.5 * (layers[i - 1] + (layers[i] if i < count else 0.0)) / dt
               for i in range(1, count + 1)]
    lower = [-conductance[i - 1] for i in range(1, count + 1)]
    upper = [-conductance[i] if i < count else 0.0 for i in range(1, count + 1)]
    diagonal = [storage[k] - lower[k] - upper[k] for k in range(count)]
    rhs = [storage[k] * float(old[k + 1]) for k in range(count)]
    rhs[0] += conductance[0] * float(surface)
    for k in range(1, count):
        multiplier = lower[k] / diagonal[k - 1]
        diagonal[k] -= multiplier * upper[k - 1]
        rhs[k] -= multiplier * rhs[k - 1]
    result = [0.0] * count
    result[-1] = rhs[-1] / diagonal[-1]
    for k in range(count - 2, -1, -1):
        result[k] = (rhs[k] - upper[k] * result[k + 1]) / diagonal[k]
    return [float(surface)] + result


def _column_start(state, profile_key, temp, nodes):
    old = [float(value) for value in state.get(profile_key, ())] or [float(temp)]
    # Missing deeper nodes take the deepest known temperature.
    return (old + [old[-1]] * nodes)[:nodes]


def conduction_profile(settings, state, temp, profile_key, urban=False, envelope=False):
    """Predict a profile without advancing the persistent previous-time state.

    Returns (profile, surface layer thickness, conductivity).
    """
    conductivity, capacity = thermal_properties(settings, urban=urban)
    dt = time_step(settings, state)
    if envelope:
        layers = max(int(settings.get("urban_wall_layers", 3)), 1)
        dz = max(float(settings.get("urban_wall_thickness", 0.5)), 1.0e-3) / layers
        bottom = float(settings.get("urban_internal_temperature_c", 24.0)) + 273.15
        old = list(state.get(profile_key, (bottom,) * (layers + 1)))
        old = (old[:layers + 1] + [bottom] * (layers + 1))[:layers + 1]
        return implicit_profile(old, temp, bottom, conductivity, capacity, dz, dt), dz, conductivity
    layers = column_layers()
    old = _column_start(state, profile_key, temp, len(layers) + 1)
    return column_profile(old, temp, conductivity, capacity, layers, dt), layers[0], conductivity


def column_surface_flux(settings, state, temp, profile_key, urban=False):
    """Surface flux into the column: top-layer conduction plus half-layer storage."""
    profile, dz, conductivity = conduction_profile(settings, state, temp, profile_key, urban=urban)
    _, capacity = thermal_properties(settings, urban=urban)
    previous = _column_start(state, profile_key, temp, len(profile))[0]
    storage = capacity * 0.5 * dz * (float(temp) - previous) / time_step(settings, state)
    return conductivity * (float(temp) - profile[1]) / dz + storage


def history_intervals(state, count, current_dt):
    """Newest-first durations between stored temperature samples.

    Missing durations in legacy hand-built states use the current interval;
    normal pipeline states explicitly carry every interval from initialization.
    """
    stored = list(state.get("thermal_history_dt_seconds", ()))
    return [max(float(stored[index] if index < len(stored) else current_dt), 1.0)
            for index in range(count)]


def surface_memory_flux(temp, history, intervals, current_dt, inertia):
    """Integrate piecewise-linear temperature changes over the last 48 h.

    Each contribution integrates (dT/dt)/sqrt(age) over its true interval.
    A partly overlapping oldest interval retains its original temperature slope.
    """
    newer = float(temp)
    elapsed = 0.0
    weighted = 0.0
    durations = [current_dt] + list(intervals)
    for previous, duration in zip(history, durations):
        end = min(elapsed + duration, HISTORY_SECONDS)
        weighted += ((newer - float(previous)) / duration
                     * (math.sqrt(end) - math.sqrt(elapsed)))
        elapsed += duration
        if elapsed >= HISTORY_SECONDS:
            break
        newer = float(previous)
    return 2.0 * inertia / math.sqrt(math.pi) * weighted
