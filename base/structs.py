from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, Mapping, MutableMapping, Optional


def _number(values: Mapping[str, Any], key: str, default: float = 0.0) -> float:
    value = values.get(key, default)
    return float(default if value in (None, "") else value)


def _numbers(values: Mapping[str, Any], key: str) -> tuple[float, ...]:
    value = values.get(key, "")
    if value in (None, ""):
        return ()
    return tuple(float(item.strip()) for item in str(value).replace(";", ",").split(",") if item.strip())


def _time_numbers(values: Mapping[str, Any], key: str) -> tuple[float, ...]:
    """Read absolute day-of-year times; every item must be a finite number.

    Wildcards such as ``*.5`` are not supported: use observation_time_mode=hour with
    observation_hours for clock hours repeated every day.
    """
    value = values.get(key, "")
    if value in (None, ""):
        return ()
    numbers = []
    for item in str(value).replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        try:
            number = float(item)
        except ValueError:
            number = float("nan")
        if number != number or number in (float("inf"), float("-inf")):
            hint = " (wildcards are not supported; use observation_time_mode=hour with observation_hours)" if "*" in item else ""
            raise ValueError(f"{key} must contain day-of-year times separated by ';'; invalid item {item!r}{hint}")
        numbers.append(number)
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
    terrain_type: int = 2
    terrain_slope: float = 0.0
    terrain_aspect: float = 0.0
    terrain_radius: float = 10.0
    terrain_height: float = 20.0
    terrain_density: float = 0.001
    terrain_base_model: str = "hom"
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
    latitude: float = 38.8555
    longitude: float = 100.3722
    dem: float = 1556.0
    slope: float = 0.0
    aspect: float = 0.0
    vza: float = 0.0
    vaa: float = 0.0
    wavelength: float = 10.5
    max_iter: int = 50
    soil_heat_method: str = "surface"
    urban_heat_method: str = "depth"
    temperature_method: str = "robust"
    stomatal_method: str = "farquhar"
    stomatal_g0: float = 0.01
    soil_heat_fraction: float = 0.25
    urban_heat_fraction: float = 0.4
    soil_temperature_method: str = "robust"
    urban_temperature_method: str = "robust"
    urban_rss: float = 3662.8562
    urban_cs: float = 1000000.0
    urban_rhos: float = 1800.0
    urban_lambda: float = 1.5
    urban_rbs: float = 10.0
    urban_temperature_c: float = 20.0
    urban_emissivity_roof: float = 0.95
    urban_emissivity_wall: float = 0.92
    urban_emissivity_street: float = 0.955
    lai: float = 0.0
    height: float = 3.0
    gleaf: float = 0.5
    lidfa: float = -0.35
    lidfb: float = -0.15
    hspot: float = 0.2
    leaf_width: float = 0.2
    row_width: float = 0.8
    row_blank: float = 0.8
    row_height: float = 1.5
    row_azimuth: float = 0.0
    crown_rad_a: float = 3.0
    crown_rad_b: float = 3.0
    crown_height: float = 3.0
    crown_density: float = 0.02
    building_density: float = 0.0015
    building_length: float = 10.0
    building_width: float = 10.0
    building_height: float = 20.0
    observation_wavelengths: tuple[float, ...] = ()
    observation_times: tuple[float, ...] = ()
    observation_time_mode: str = "all"
    observation_vza: tuple[float, ...] = ()
    observation_vaa: tuple[float, ...] = ()
    observation_direction_mode: str = "file"
    directions_file: Optional[Path] = None
    observation_spectral_mode: str = "file"
    spectrals_file: Optional[Path] = None
    output_intermediate: bool = False
    raw: Dict[str, Any] = field(default_factory=dict, repr=False)

    @classmethod
    def from_config(cls, config: Mapping[str, Any], base: Path) -> "InputStructure":
        def path(key: str, default: str) -> Path:
            value = Path(str(config.get(key, default)))
            return value if value.is_absolute() else base / value

        model = str(config.get("surface_model", "plane")).strip().lower().replace("-", "_")
        if model == "slope":
            model = "terrain"
        vegetation = str(config.get("vegetation_model", config.get("vegetation_type", "hom"))).strip().lower()
        def numbers(*keys: str) -> tuple[float, ...]:
            for key in keys:
                result = _numbers(config, key)
                if result:
                    return result
            return ()
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
        direction_source = config.get("direction_source", config.get("observation_direction_mode", "file"))
        spectral_source = config.get("spectral_source", config.get("observation_spectral_mode", "direct"))
        return cls(
            surface_model=model,
            vegetation_model=vegetation,
            terrain_type=1 if str(config.get("surface_model", "")).strip().lower() == "slope" else int(float(config.get("terrain_type", 2))),
            terrain_slope=_number(config, "terrain_slope"),
            terrain_aspect=_number(config, "terrain_aspect"),
            terrain_radius=_number(config, "terrain_radius", 10.0),
            terrain_height=_number(config, "terrain_height", 20.0),
            terrain_density=_number(config, "terrain_density", 0.001),
            terrain_base_model=str(config.get("terrain_base_model", "hom")),
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
            latitude=_number(config, "latitude", 38.8555),
            longitude=_number(config, "longitude", 100.3722),
            dem=_number(config, "dem", 1556.0),
            slope=_number(config, "slope"),
            aspect=_number(config, "aspect", _number(config, "slope_aspect")),
            vza=(numbers("vza", "observation_vza") or (0.0,))[0],
            vaa=(numbers("vaa", "observation_vaa") or (0.0,))[0],
            wavelength=_number(config, "wavelength", 10.5),
            max_iter=int(float(config.get("max_iter", 50))),
            soil_heat_method=str(config.get("soil_heat_method", "surface")).strip().lower(),
            urban_heat_method=str(config.get("urban_heat_method", "depth")).strip().lower(),
            temperature_method=str(config.get("temperature_method", "robust")).strip().lower(),
            stomatal_method=str(config.get("stomatal_method", "farquhar")).strip().lower(),
            stomatal_g0=_number(config, "stomatal_g0", 0.01),
            soil_heat_fraction=_number(config, "soil_heat_fraction", 0.25),
            urban_heat_fraction=_number(config, "urban_heat_fraction", 0.4),
            soil_temperature_method=str(config.get("soil_temperature_method") or config.get("temperature_method") or "robust").strip().lower(),
            urban_temperature_method=str(config.get("urban_temperature_method") or config.get("temperature_method") or "robust").strip().lower(),
            urban_rss=_number(config, "urban_rss", 3662.8562),
            urban_cs=_number(config, "urban_cs", 1000000.0),
            urban_rhos=_number(config, "urban_rhos", 1800.0),
            urban_lambda=_number(config, "urban_lambda", 1.5),
            urban_rbs=_number(config, "urban_rbs", 10.0),
            urban_temperature_c=_number(config, "urban_temperature_c", 20.0),
            urban_emissivity_roof=_number(config, "urban_emissivity_roof", 0.95),
            urban_emissivity_wall=_number(config, "urban_emissivity_wall", 0.92),
            urban_emissivity_street=_number(config, "urban_emissivity_street", 0.955),
            output_intermediate=_boolean(config, "output_intermediate", False),
            lai=_number(config, "lai"),
            height=_number(config, "height", 3.0),
            gleaf=_number(config, "gleaf", 0.5),
            lidfa=_number(config, "lidfa", -0.35),
            lidfb=_number(config, "lidfb", -0.15),
            hspot=_number(config, "hspot", 0.2),
            leaf_width=_number(config, "leaf_width", 0.2),
            row_width=_number(config, "row_width", _number(config, "ridge_width", 0.8)),
            row_blank=_number(config, "row_blank", _number(config, "ridge_blank", 0.8)),
            row_height=_number(config, "row_height", _number(config, "ridge_height", 1.5)),
            row_azimuth=_number(config, "row_azimuth", _number(config, "ridge_azimuth")),
            crown_rad_a=_number(config, "crown_rad_a", 3.0),
            crown_rad_b=_number(config, "crown_rad_b", 3.0),
            crown_height=_number(config, "crown_height", 3.0),
            crown_density=_number(config, "crown_density", 0.02),
            building_density=_number(config, "building_density", 0.0015),
            building_length=_number(config, "building_length", 10.0),
            building_width=_number(config, "building_width", 10.0),
            building_height=_number(config, "building_height", 20.0),
            observation_wavelengths=numbers("wavelengths", "observation_wavelengths"),
            observation_times=_time_numbers(config, "observation_times"),
            observation_time_mode=str(config.get("observation_time_mode", "all")).strip().lower(),
            observation_vza=numbers("vza", "observation_vza"),
            observation_vaa=numbers("vaa", "observation_vaa"),
            observation_direction_mode=str(direction_source).strip().lower(),
            directions_file=path("directions_file", "data/directions/observation_directions.txt"),
            observation_spectral_mode=str(spectral_source).strip().lower(),
            spectrals_file=path("spectrals_file", "data/spectrals/observation_spectrals_tirteb.txt"),
            raw=dict(config),
        )

    def as_config(self) -> Dict[str, Any]:
        return dict(self.raw)

    def active_structure(self) -> Dict[str, Any]:
        common = {"lai": self.lai, "height": self.height, "leaf_width": self.leaf_width, "hspot": self.hspot}
        if self.vegetation_model == "bare":
            return {"lai": 0.0}
        if self.vegetation_model == "row":
            return {**common, "row_width": self.row_width, "row_blank": self.row_blank,
                    "row_height": self.row_height, "row_azimuth": self.row_azimuth}
        if self.vegetation_model == "crown":
            return {**common, "crown_rad_a": self.crown_rad_a, "crown_rad_b": self.crown_rad_b,
                    "crown_height": self.crown_height, "crown_density": self.crown_density,
                    "gleaf": self.gleaf, "lidfa": self.lidfa, "lidfb": self.lidfb}
        if self.surface_model == "terrain":
            return {**common, "terrain_type": self.terrain_type,
                    "terrain_slope": self.terrain_slope, "terrain_aspect": self.terrain_aspect,
                    "terrain_radius": self.terrain_radius, "terrain_height": self.terrain_height,
                    "terrain_density": self.terrain_density}
        if self.surface_model == "urban":
            return {"building_density": self.building_density,
                    "building_length": self.building_length, "building_width": self.building_width,
                    "building_height": self.building_height}
        return common


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
            "canopy": ("lai", "stand_density", "height", "width", "gleaf", "lidfa", "lidfb", "hspot", "leaf_width"),
            "leaf": ("leaf_reflectance", "leaf_transmittance", "leaf_emissivity"),
            "biochemistry": ("vcmax", "ballberry_m", "leaf_type", "kv", "rdparam", "tparam", "tyear", "beta", "knpqs", "qls", "stressfactor", "stomatal_method", "stomatal_g0"),
            "soil_optical": ("soil_reflectance", "soil_emissivity", "hapke_enable", "hapke_K"),
            "soil_thermal": ("soil_rss", "soil_cs", "soil_rhos", "soil_lambda", "soil_rbs", "soil_moisture", "soil_temperature_c", "satwater", "soil_tir_water_effect", "soil_tir_emissivity_wet", "soil_tir_moisture_dry", "soil_tir_moisture_wet"),
            "urban_thermal": ("urban_rss", "urban_cs", "urban_rhos", "urban_lambda", "urban_rbs", "urban_temperature_c", "urban_wall_thickness", "urban_wall_layers", "urban_internal_temperature_c", "urban_wall_diffuse_factor", "urban_emissivity_roof", "urban_emissivity_wall", "urban_emissivity_street"),
            "urban_spectrum": ("urban_spectrum_source", "urban_spectrum_file", "urban_reflectance_roof_spectrum", "urban_reflectance_wall_spectrum", "urban_reflectance_street_spectrum", "urban_emissivity_roof_spectrum", "urban_emissivity_wall_spectrum", "urban_emissivity_street_spectrum"),
            "run": ("soil_heat_method", "urban_heat_method", "soil_temperature_method", "urban_temperature_method", "soil_heat_fraction", "urban_heat_fraction", "temperature_method"),
            "aerodynamic": ("aero_zo", "aero_d", "aero_rbc", "aero_cr", "aero_cd", "aero_cd1", "aero_psicor", "aero_cssoil", "aero_rwc", "aero_rbs"),
            "atmosphere": ("z", "ea", "ca", "oa", "shortwave_diffuse_fraction", "sky_temperature_k"),
            "spectra": ("wavelengths", "fesun", "fesky", "leaf_reflectance_spectrum", "leaf_transmittance_spectrum", "soil_reflectance_spectrum"),
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
