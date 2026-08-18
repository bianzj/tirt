"""Directional soil emissivity helpers based on a Hapke/Thapke form."""

import numpy as np


def hapke_emissivity(view_zenith, K=1.0, albedo=0.05):
    """Return directional thermal emissivity for a soil albedo."""
    angle = np.asarray(view_zenith, dtype=float)
    coefficient = max(float(K), 1.0e-8)
    albedo = np.clip(np.asarray(albedo, dtype=float), 0.0, 0.999999)
    base_emissivity = np.sqrt(1.0 - albedo)
    term = 2.0 * np.cos(np.deg2rad(angle)) / coefficient
    denominator = 1.0 + term * base_emissivity
    result = base_emissivity * (1.0 + term) / np.maximum(denominator, 1.0e-12)
    return np.clip(result, 0.0, 1.0)


def hapke_thermal_emissivity(view_zenith, K=1.0, albedo=0.05):
    """Convenience wrapper for the thermal Hapke interface."""
    return hapke_emissivity(view_zenith, K, albedo)


__all__ = ["hapke_emissivity", "hapke_thermal_emissivity"]
