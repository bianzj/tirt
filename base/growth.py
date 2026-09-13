def update_canopy(settings, lai=None):
    if lai is not None:
        settings["lai"] = max(float(lai), 0.0)
    return settings
