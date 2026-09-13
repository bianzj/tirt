import math


def clamp(value, lower, upper):
    return max(lower, min(upper, value))


def finite_or(value, fallback):
    return value if math.isfinite(value) else fallback
