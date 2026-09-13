import importlib
import sys
import types


def _optional_dependencies():
    try:
        importlib.import_module("scipy.integrate")
    except ImportError:
        scipy = types.ModuleType("scipy")
        integrate = types.ModuleType("scipy.integrate")
        scipy.integrate = integrate
        sys.modules["scipy"] = scipy
        sys.modules["scipy.integrate"] = integrate
    try:
        importlib.import_module("matplotlib.pyplot")
    except ImportError:
        matplotlib = types.ModuleType("matplotlib")
        pyplot = types.ModuleType("matplotlib.pyplot")
        matplotlib.pyplot = pyplot
        sys.modules["matplotlib"] = matplotlib
        sys.modules["matplotlib.pyplot"] = pyplot
    try:
        importlib.import_module("mpl_toolkits.mplot3d")
    except ImportError:
        mpl = types.ModuleType("mpl_toolkits")
        mplot3d = types.ModuleType("mpl_toolkits.mplot3d")
        mplot3d.Axes3D = object
        mpl.mplot3d = mplot3d
        sys.modules["mpl_toolkits"] = mpl
        sys.modules["mpl_toolkits.mplot3d"] = mplot3d


_optional_dependencies()

import numpy as np

if not hasattr(np, "float_"):
    np.float_ = np.float64
if not hasattr(np, "int"):
    np.int = int

from .pipeline import run_node
from .structs import (
    DefinedStructure, FluxStructure, InputStructure, ObservationStructure,
    RTStructure,
)

__all__ = [
    "run_node", "InputStructure", "DefinedStructure", "RTStructure",
    "FluxStructure", "ObservationStructure",
]
