"""Two-scale forest geometry applied to thermal emission (wavelength in um).

clumping_dual uses whole-crown analytic averages; volume uses the inherited
FRT crown cubature. Only geometric probabilities are shared with optics.
The existing analytic thermal multiple-scattering closure is retained,
with its contribution limited to the available directional emissivity budget.
"""
import numpy as np

from .crown_parametric import crown_fractions
from .emissivity import canopy_scatter
from .optical_crown import fractions
from .physicsF import inv_planck, planck


MODES = ("legacy", "clumping_dual", "volume", "parametric")


def geometry_mode(settings):
    mode = str(settings.get("thermal_crown_geometry", "legacy")).strip().lower()
    if mode not in MODES:
        raise ValueError(f"thermal_crown_geometry must be one of {MODES}.")
    return mode


def is_forest(settings):
    model = str(settings.get("vegetation_model", "")).strip().lower()
    return model in ("crown", "forest") or int(settings.get("surface_type", 0)) in (1, 2, 3, 4, 5)


def component_fractions(settings, sza, vza, raa, sunlit=True):
    """Return (4, N) soil-sun/shade, leaf-sun/shade probabilities.

    height/width are ellipsoid vertical/horizontal semiaxes. The volume
    adapter uses crown bottom zero, matching the existing optical adapter.
    thermal_hotspot_inner/outer are supported only by clumping_dual.
    """
    mode = geometry_mode(settings)
    if mode == "legacy":
        raise ValueError("Use the original thermal entry point for legacy geometry.")
    shape = tuple(float(settings[k]) for k in ("lai", "stand_density", "height", "width", "hspot"))
    lai, density, height, width, hspot = shape
    if not np.isfinite(shape).all() or min(lai, density, hspot) < 0 or min(height, width) <= 0:
        raise ValueError("Finite nonnegative LAI/density/hspot and positive crown semiaxes required.")
    scales = tuple(float(settings.get(key) if settings.get(key) not in (None, "") else default)
                   for key, default in (("thermal_hotspot_inner", hspot),
                                        ("thermal_hotspot_outer", width / height)))
    if not np.isfinite(scales).all() or min(scales) < 0:
        raise ValueError("Thermal hotspot widths must be finite and nonnegative.")
    if mode == "volume" and any(settings.get(k) not in (None, "") for k in
                                ("thermal_hotspot_inner", "thermal_hotspot_outer")):
        raise ValueError("thermal_hotspot_inner/outer require clumping_dual; volume uses hspot.")
    angles = np.broadcast_arrays(np.atleast_1d(sza), np.atleast_1d(vza), np.atleast_1d(raa))
    if not all(np.isfinite(a).all() for a in angles):
        raise ValueError("Thermal angles must be finite.")
    solar, view, azimuth = (np.asarray(a, dtype=float).ravel() for a in angles)
    solar, view = np.clip(solar, 0, 89.9), np.clip(np.abs(view), 0, 89.9)
    azimuth = np.abs((azimuth + 180) % 360 - 180)
    if lai == 0 or density == 0:
        values = np.zeros((4, solar.size))
        values[0] = 1
    elif mode == "parametric":
        base = settings.get("crown_base")
        base = height if base in (None, "") else float(base)
        values = crown_fractions(lai, density, hspot, height, width, base, view, solar, azimuth)
    else:
        values = np.asarray([fractions(*shape, float(s), float(v), float(a), mode,
                                       scales if mode == "clumping_dual" else None)
                             for s, v, a in zip(solar, view, azimuth)]).T
    if not sunlit:
        values[1] += values[0]
        values[3] += values[2]
        values[[0, 2]] = 0
    return values


def thermal_emission(settings, temperatures, wavelength, sza, vza, raa, sunlit=True):
    """Return radiance, brightness temperature and four geometric fractions.

    L = sum_i [epsilon_i * p_i + scatter_i] * B_lambda(T_i).
    No optical BRF correction is applied. Temperature contrasts, rather than
    solar illumination itself, determine the sign of a thermal hotspot.
    """
    p = component_fractions(settings, sza, vza, raa, sunlit)
    es, el = float(settings["soil_emissivity"]), float(settings["leaf_emissivity"])
    t = np.asarray(temperatures, dtype=float).reshape(4)
    if not np.isfinite(t).all() or np.any(t <= 0) or not np.isfinite(wavelength) or wavelength <= 0:
        raise ValueError("Temperatures (K) and wavelength (um) must be finite and positive.")
    if not all(np.isfinite(e) and 0 <= e <= 1 for e in (es, el)):
        raise ValueError("Thermal emissivities must be in [0, 1].")
    scatter = np.zeros_like(p)
    if settings["lai"] > 0 and settings["stand_density"] > 0:
        s, v, _ = np.broadcast_arrays(np.atleast_1d(sza), np.atleast_1d(vza), np.atleast_1d(raa))
        scatter = canopy_scatter(
            "crown", float(settings["lai"]), np.array([el]), np.array([es]),
            np.clip(np.abs(v.astype(float).ravel()), 0, 89.9),
            np.clip(s.astype(float).ravel(), 0, 89.9),
            float(settings["stand_density"]), float(settings["height"]), float(settings["width"]))
        scatter = np.asarray(scatter).reshape(4, -1)
        if not sunlit:
            scatter[1] += scatter[0]
            scatter[3] += scatter[2]
            scatter[[0, 2]] = 0
    direct = p * np.array([es, es, el, el])[:, None]
    # The inherited homogeneous-canopy scatter closure can exceed the
    # remaining budget for sparse crowns over a highly emissive substrate.
    # Limit the scatter only; preserve its sun/shade ratio and all geometry.
    # This guard is not a replacement for a forest photon-transport closure.
    budget = np.maximum(0., 1. - direct.sum(axis=0))
    total_scatter = scatter.sum(axis=0)
    scale = np.ones_like(budget)
    np.divide(budget, total_scatter, out=scale, where=total_scatter > 0)
    scatter *= np.minimum(1., scale)[None, :]
    weights = direct + scatter
    radiance = np.sum(weights * planck(wavelength, t)[:, None], axis=0)
    return radiance, inv_planck(wavelength, radiance), p
