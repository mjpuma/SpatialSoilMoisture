"""Spatial soil moisture and evapotranspiration.

The model is the last part of Puma, M.J. (2006), Ph.D. dissertation,
Princeton University, and of the public text of that chapter,
arXiv:1606.05256.

Start with `simulate_pair` for the uniform-rain, grass-and-woody case,
or with `daily_step` for one root zone. Depths are millimeters and time
is days. `tests/` is the check that the drying solution, the rainfall
moments, the woody cover, and the water balance are doing what the
dissertation specifies. Run them from the repository root:

    PYTHONPATH=src python -m pytest
"""

from .parameters import (
    GRASS,
    RAIN,
    SOIL,
    WOODY,
    PlantParams,
    RainParams,
    SoilParams,
    effective_plant,
    mean_storm_depth,
)
from .point import apply_rain, daily_step, dry, et_rate, leakage_rate

__all__ = [
    "GRASS",
    "RAIN",
    "SOIL",
    "WOODY",
    "PlantParams",
    "RainParams",
    "SoilParams",
    "apply_rain",
    "daily_step",
    "dry",
    "effective_plant",
    "et_rate",
    "leakage_rate",
    "mean_storm_depth",
]
