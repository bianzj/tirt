from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, Mapping, MutableMapping, Optional


def _number(values: Mapping[str, Any], key: str, default: float = 0.0) -> float:
    value = values.get(key, default)
    return float(default if value in (None, "") else value)


def _time_numbers(values: Mapping[str, Any], key: str) -> tuple[float, ...]:
    """Read numeric times while allowing wildcard hour syntax such as *.5."""
    value = values.get(key, "")
    if value in (None, ""):
        return ()
    numbers = []
    for item in str(value).replace(";", ",").split(","):
        item = item.strip()
        if not item or "*" in item:
            continue
        try:
            numbers.append(float(item))
        except ValueError:
            continue
    return tuple(numbers)


def _boolean(values: Mapping[str, Any], key: str, default: bool = False) -> bool:
    value = values.get(key, default)
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off"}:
        return False
    return default


class MappingStruct(MutableMapping[str, Any]):
    """Dataclass container that remains compatible with the legacy dict API."""

    values: Dict[str, Any]

    def __getitem__(self, key: str) -> Any:
        return self.values[key]

    def __setitem__(self, key: str, value: Any) -> None:
        self.values[key] = value

    def __delitem__(self, key: str) -> None:
        del self.values[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self.values)

    def __len__(self) -> int:
        return len(self.values)

    def get(self, key: str, default: Any = None) -> Any:
        return self.values.get(key, default)

    def update(self, *args: Any, **kwargs: Any) -> None:
        self.values.update(*args, **kwargs)

    def as_dict(self) -> Dict[str, Any]:
        return dict(self.values)


@dataclass
class Angle:
    vza: float = 0.0
    vaa: float = 0.0
    sza: float = 0.0
    saa: float = 0.0


@dataclass
class Thermal:
    Tss: float = 273.17
    Tsh: float = 273.17
    Tls: float = 273.17
    Tlh: float = 273.17


@dataclass
class Resistance:
    L: float = 1e9
    ustar: float = 0.0
    raa_leaf: float = 0.0
    raa_soil: float = 0.0
    raa_roof: float = 0.0
    raa_wall: float = 0.0
    raa_street: float = 0.0


@dataclass
class InputStructure:
    """All values read from input.csv, including active surface structure."""

    surface_model: str = "plane"
    vegetation_model: str = "hom"
    defined_path: Optional[Path] = None
    meteo_path: Optional[Path] = None
    cases_dir: Optional[Path] = None
    case_name: str = "homogeneous"
    case_dir: Optional[Path] = None
    output_path: Optional[Path] = None
    year: int = 2024
    doy: int = 172
    utc_offset_hours: float = 0.0
    start_time: float = float("-inf")
    end_time: float = float("inf")
    latitude: float = 39.9
    longitude: float = 116.4
    dem: float = 0.0
    max_iter: int = 100
    soil_heat_method: str = "surface"
    urban_heat_method: str = "depth"
    temperature_method: str = "robust"
    stomatal_method: str = "farquhar"
    stomatal_g0: float = 0.01
    soil_heat_fraction: float = 0.25
    urban_heat_fraction: float = 0.4
    soil_depth_step: float = 0.05
    soil_temperature_method: str = "robust"
    urban_temperature_method: str = "robust"
    urban_rss: float = 3662.8562
    urban_cs: float = 1000000.0
    urban_rhos: float = 1800.0
    urban_lambda: float = 1.5
    urban_rbs: float = 10.0
    urban_temperature_c: float = 20.0
    urban_depth_step: float = 0.05
    observation_times: tuple[float, ...] = ()
    observation_time_mode: str = "all"
    output_intermediate: bool = False
    raw: Dict[str, Any] = field(default_factory=dict, repr=False)

    @classmethod
    def from_config(cls, config: Mapping[str, Any], base: Path) -> "InputStructure":
        def path(key: str, default: str) -> Path:
            value = Path(str(config.get(key, default)))
            return value if value.is_absolute() else base / value

        model = str(config.get("surface_model", "plane")).strip().lower().replace("-", "_")
        vegetation = str(config.get("vegetation_model", config.get("vegetation_type", "hom"))).strip().lower()
        # A normal run is always stored below the configured cases directory.
        # case_result_dir remains an internal override used by --compare.
        if config.get("case_result_dir") not in (None, ""):
            result_dir = config["case_result_dir"]
        else:
            cases_root = path("cases_dir", "cases")
            case_name = str(config.get("case_name", "homogeneous")).strip()
            if not case_name:
                case_name = "homogeneous"
            result_dir = cases_root / case_name
        return cls(
            surface_model=model,
            vegetation_model=vegetation,
            defined_path=None,
            meteo_path=path("meteo", "data/defined/meteo.txt"),
            cases_dir=path("cases_dir", "cases"),
            case_name=str(config.get("case_name", "homogeneous")),
            case_dir=Path(str(result_dir)) if Path(str(result_dir)).is_absolute() else base / str(result_dir),
            output_path=None,
            year=int(float(config.get("year", 2024))),
            doy=int(float(config.get("doy", 172))),
            utc_offset_hours=_number(config, "utc_offset_hours", _number(config, "utc_offset", 0.0)),
            start_time=_number(config, "start_time", float("-inf")),
            end_time=_number(config, "end_time", float("inf")),
            latitude=_number(config, "latitude", 39.9),
            longitude=_number(config, "longitude", 116.4),
            dem=_number(config, "dem"),
            max_iter=int(float(config.get("max_iter", 100))),
            soil_heat_method=str(config.get("soil_heat_method", "surface")).strip().lower(),
            urban_heat_method=str(config.get("urban_heat_method", "depth")).strip().lower(),
            temperature_method=str(config.get("temperature_method", "robust")).strip().lower(),
            stomatal_method=str(config.get("stomatal_method", "farquhar")).strip().lower(),
            stomatal_g0=_number(config, "stomatal_g0", 0.01),
            soil_heat_fraction=_number(config, "soil_heat_fraction", 0.25),
            urban_heat_fraction=_number(config, "urban_heat_fraction", 0.4),
            soil_depth_step=_number(config, "soil_depth_step", 0.05),
            soil_temperature_method=str(config.get("soil_temperature_method", "robust")).strip().lower(),
            urban_temperature_method=str(config.get("temperature_method", config.get("urban_temperature_method", "robust"))).strip().lower(),
            urban_rss=_number(config, "urban_rss", 3662.8562),
            urban_cs=_number(config, "urban_cs", 1000000.0),
            urban_rhos=_number(config, "urban_rhos", 1800.0),
            urban_lambda=_number(config, "urban_lambda", 1.5),
            urban_rbs=_number(config, "urban_rbs", 10.0),
            urban_temperature_c=_number(config, "urban_temperature_c", 20.0),
            urban_depth_step=_number(config, "urban_depth_step", 0.05),
            output_intermediate=_boolean(config, "output_intermediate", False),
            observation_times=_time_numbers(config, "observation_times"),
            observation_time_mode=str(config.get("observation_time_mode", "all")).strip().lower(),
            raw=dict(config),
        )

    def as_config(self) -> Dict[str, Any]:
        return dict(self.raw)


@dataclass
class DefinedStructure(MappingStruct):
    """Defined-file parameters/constants with legacy mapping compatibility."""

    values: Dict[str, Any] = field(default_factory=dict)
    canopy: Dict[str, Any] = field(default_factory=dict)
    leaf: Dict[str, Any] = field(default_factory=dict)
    biochemistry: Dict[str, Any] = field(default_factory=dict)
    soil_optical: Dict[str, Any] = field(default_factory=dict)
    soil_thermal: Dict[str, Any] = field(default_factory=dict)
    urban_thermal: Dict[str, Any] = field(default_factory=dict)
    urban_spectrum: Dict[str, Any] = field(default_factory=dict)
    run: Dict[str, Any] = field(default_factory=dict)
    aerodynamic: Dict[str, Any] = field(default_factory=dict)
    atmosphere: Dict[str, Any] = field(default_factory=dict)
    spectra: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> "DefinedStructure":
        data = dict(values)
        groups = {
            "canopy": ("lai", "canopy_height", "gleaf", "lidfa", "lidfb", "hspot", "leaf_width", "crown_rad_a", "crown_rad_b", "crown_height", "crown_density", "row_width", "row_blank", "row_height", "row_azimuth"),
            "leaf": ("leaf_reflectance", "leaf_transmittance", "leaf_emissivity"),
            "biochemistry": ("vcmax", "ballberry_m", "leaf_type", "kv", "rdparam", "tparam", "tyear", "beta", "knpqs", "qls", "stressfactor", "stomatal_method", "stomatal_g0"),
            "soil_optical": ("soil_reflectance", "soil_emissivity", "hapke_enable", "hapke_K"),
            "soil_thermal": ("soil_rss", "soil_cs", "soil_rhos", "soil_lambda", "soil_rbs", "soil_moisture", "soil_temperature_c", "satwater", "soil_depth_step"),
            "urban_thermal": ("urban_rss", "urban_cs", "urban_rhos", "urban_lambda", "urban_rbs", "urban_temperature_c", "urban_depth_step", "urban_wall_thickness", "urban_wall_layers", "urban_internal_temperature_c", "urban_wall_diffuse_factor", "urban_emissivity_roof", "urban_emissivity_wall", "urban_emissivity_street"),
            "urban_spectrum": ("urban_reflectance_roof_spectrum", "urban_reflectance_wall_spectrum", "urban_reflectance_street_spectrum"),
            "run": ("soil_heat_method", "urban_heat_method", "soil_temperature_method", "urban_temperature_method", "soil_heat_fraction", "urban_heat_fraction", "temperature_method"),
            "aerodynamic": ("aero_zo", "aero_d", "aero_rbc", "aero_cr", "aero_cd", "aero_cd1", "aero_psicor", "aero_cssoil", "aero_rwc", "aero_rbs"),
            "atmosphere": ("z", "ca", "oa", "shortwave_diffuse_fraction"),
            "spectra": ("wavelengths", "fesun", "fesky", "leaf_reflectance_spectrum", "leaf_transmittance_spectrum", "soil_reflectance_spectrum", "tir_wavelengths", "tir_emissivity"),
        }
        return cls(data, **{name: {key: data[key] for key in keys if key in data} for name, keys in groups.items()})


@dataclass
class RTStructure(MappingStruct):
    """One wavelength/view radiative-transfer result."""

    values: Dict[str, Any] = field(default_factory=dict)
    wavelength_um: float = float("nan")
    reflectance: float = float("nan")
    transmittance: float = float("nan")
    emissivity: float = float("nan")
    radiance: float = float("nan")
    brightness_temperature_K: float = float("nan")

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> "RTStructure":
        data = dict(values)
        return cls(data, float(data.get("wavelength_um", float("nan"))),
                   float(data.get("reflectance", float("nan"))),
                   float(data.get("transmittance", float("nan"))),
                   float(data.get("emissivity", float("nan"))),
                   float(data.get("radiance", float("nan"))),
                   float(data.get("brightness_temperature_K", float("nan"))))


@dataclass
class FluxStructure(MappingStruct):
    """Energy-balance state and component intermediate variables."""

    values: Dict[str, Any] = field(default_factory=dict)
    thermal: Thermal = field(default_factory=Thermal)
    resistance: Resistance = field(default_factory=Resistance)

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> "FluxStructure":
        data = dict(values)
        thermal = Thermal(*(float(data.get(key, default)) for key, default in (
            ("Tss", 273.17), ("Tsh", 273.17), ("Tls", 273.17), ("Tlh", 273.17))))
        resistance = Resistance(*(float(data.get(key, default)) for key, default in (
            ("L", 1e9), ("ustar", 0.0), ("raa_leaf", 0.0), ("raa_soil", 0.0),
            ("raa_roof", 0.0), ("raa_wall", 0.0), ("raa_street", 0.0))))
        return cls(data, thermal, resistance)

    def sync_structures(self) -> None:
        self.thermal = Thermal(self.get("Tss", 273.17), self.get("Tsh", 273.17), self.get("Tls", 273.17), self.get("Tlh", 273.17))
        self.resistance = Resistance(*(self.get(key, default) for key, default in (
            ("L", 1e9), ("ustar", 0.0), ("raa_leaf", 0.0), ("raa_soil", 0.0),
            ("raa_roof", 0.0), ("raa_wall", 0.0), ("raa_street", 0.0))))


@dataclass
class ObservationStructure:
    """One selected time/view observation with wavelength columns."""

    node: int
    time: float
    model: str
    surface_model: str
    vza: float
    vaa: float
    sza: float
    saa: float
    results: Dict[str, Any] = field(default_factory=dict)

    def as_row(self) -> Dict[str, Any]:
        return {
            "node": self.node, "time": self.time, "veg": self.model,
            "surface": self.surface_model, "vza": self.vza, "vaa": self.vaa,
            "sza": self.sza, "saa": self.saa, **self.results,
        }


@dataclass
class PixelState:
    angle: Angle = field(default_factory=Angle)
    thermal: Thermal = field(default_factory=Thermal)
    resistance: Resistance = field(default_factory=Resistance)


__all__ = [
    "Angle", "Thermal", "Resistance", "InputStructure", "DefinedStructure",
    "RTStructure", "FluxStructure", "ObservationStructure", "PixelState",
]
