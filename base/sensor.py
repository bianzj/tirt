from dataclasses import dataclass


@dataclass
class Observation:
    vza: float = 0.0
    vaa: float = 0.0
    wavelength: float = 10.5
