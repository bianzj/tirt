import math

import numpy as np

from .physicsF import slope1


def _refraction(elevation):
    if elevation > 85:
        correction = 0.0
    else:
        tangent = math.tan(math.radians(elevation))
        if elevation > 5:
            correction = 58.1 / tangent - 0.07 / tangent**3 + 0.000086 / tangent**5
        elif elevation > -0.575:
            correction = 1735 + elevation * (-518.2 + elevation * (103.4 + elevation * (-12.79 + elevation * 0.711)))
        else:
            correction = -20.774 / tangent
        correction /= 3600.0
    return correction


def solar_angles(year, doy, time_value, latitude, longitude, timezone=0.0):
    """NOAA SPA-equivalent solar zenith/azimuth calculation."""
    fraction = time_value - math.floor(time_value)
    hour = fraction * 24.0
    gamma = 2.0 * math.pi / 365.0 * (doy - 1 + (hour - 12.0) / 24.0)
    equation = 229.18 * (
        0.000075 + 0.001868 * math.cos(gamma)
        - 0.032077 * math.sin(gamma) - 0.014615 * math.cos(2 * gamma)
        - 0.040849 * math.sin(2 * gamma)
    )
    declination = (
        0.006918 - 0.399912 * math.cos(gamma) + 0.070257 * math.sin(gamma)
        - 0.006758 * math.cos(2 * gamma) + 0.000907 * math.sin(2 * gamma)
        - 0.002697 * math.cos(3 * gamma) + 0.00148 * math.sin(3 * gamma)
    )
    true_solar_minutes = hour * 60 + equation + 4 * longitude - 60 * timezone
    hour_angle = true_solar_minutes / 4 - 180
    while hour_angle < -180:
        hour_angle += 360
    while hour_angle > 180:
        hour_angle -= 360

    lat = math.radians(latitude)
    dec = declination
    ha = math.radians(hour_angle)
    csz = math.sin(lat) * math.sin(dec) + math.cos(lat) * math.cos(dec) * math.cos(ha)
    csz = min(1.0, max(-1.0, csz))
    zenith = math.degrees(math.acos(csz))
    denominator = math.cos(lat) * math.sin(math.radians(zenith))
    if abs(denominator) > 0.001:
        az = (math.sin(lat) * math.cos(math.radians(zenith)) - math.sin(dec)) / denominator
        az = min(1.0, max(-1.0, az))
        azimuth = 180 - math.degrees(math.acos(az))
        if hour_angle > 0:
            azimuth = -azimuth
    else:
        azimuth = 180 if latitude > 0 else 0
    if azimuth < 0:
        azimuth += 360
    return zenith - _refraction(90 - zenith), azimuth


def terrain_solar_angles(sza, saa, slope, aspect):
    """Convert global solar angles to TIRT slope-local angles."""
    local_sza, local_saa, _ = slope_direction(sza, saa, slope, aspect)
    return local_sza, local_saa


def slope_direction(zenith, azimuth, slope, aspect):
    """Transform one direction with the same rotation used by TIRT.slope1."""
    zenith_array = np.asarray([float(zenith)], dtype=float)
    azimuth_array = np.asarray([float(azimuth)], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        local_zenith, local_azimuth = slope1(
            zenith_array, azimuth_array, float(slope), float(aspect)
        )
    rd = math.pi / 180.0
    incidence = (
        math.cos(float(slope) * rd) * math.cos(float(zenith) * rd)
        + math.sin(float(slope) * rd)
        * math.sin(float(zenith) * rd)
        * math.cos((float(azimuth) - float(aspect)) * rd)
    )
    return float(local_zenith[0]), float(local_azimuth[0]), incidence >= 0.0


__all__ = ["solar_angles", "terrain_solar_angles", "slope_direction"]
