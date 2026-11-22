"""La Copita parameters, in millimeters and days.

Import `SOIL`, `GRASS`, `WOODY`, and `RAIN` and pass them to `daily_step`
or `simulate_pair`. Relative saturation is water volume per void volume.

`effective_plant(0.5)` is the single column used when grass and woody
plants each cover half the area. Its root depth is 700 mm and its maximum
evapotranspiration is 4.59 mm/d, which is the check against Table 1 of
the dissertation chapter. `mean_storm_depth()` is about 15.3 mm. The
dissertation table rounds that to 15 mm.

The Fortran input files store the same quantities in centimeters.
`Zr = 40` there is 40 cm, which is `GRASS.zr` of 400 mm here.
`Ksat = 82.2` cm/d there is `SOIL.ksat` of 822 mm/d here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SoilParams:
    porosity: float
    ksat: float
    b: float
    s_fc: float
    s_h: float

    def __post_init__(self) -> None:
        if not 0.0 < self.s_h < self.s_fc < 1.0:
            raise ValueError("need 0 < s_h < s_fc < 1")
        if self.porosity <= 0.0 or self.ksat <= 0.0:
            raise ValueError("porosity and ksat must be positive")


@dataclass(frozen=True)
class PlantParams:
    s_w: float
    s_star: float
    zr: float
    et_max: float
    e_w: float
    interception: float
    name: str = ""

    def __post_init__(self) -> None:
        if not 0.0 < self.s_w < self.s_star < 1.0:
            raise ValueError("need 0 < s_w < s_star < 1")
        if self.zr <= 0.0 or self.et_max < 0.0 or self.e_w < 0.0:
            raise ValueError("root depth and evapotranspiration rates must be non-negative")
        if self.e_w > self.et_max:
            raise ValueError("e_w cannot exceed et_max")


@dataclass(frozen=True)
class RainParams:
    """Spatial storm model, equations (4)–(7).

    ``cell_density`` is cells per square meter. ``cell_scale`` is a' in meters.
    """

    storm_rate: float
    mean_cell_depth: float
    cell_density: float
    cell_scale: float


# Table 3. Sandy loam. Ksat is 822 mm/d.
SOIL = SoilParams(porosity=0.43, ksat=822.0, b=4.9, s_fc=0.56, s_h=0.14)

# Table 1. Paspalum setaceum.
GRASS = PlantParams(
    s_w=0.167,
    s_star=0.370,
    zr=400.0,
    et_max=4.76,
    e_w=0.13,
    interception=1.0,
    name="Paspalum setaceum",
)

# Table 1. Prosopis glandulosa.
WOODY = PlantParams(
    s_w=0.180,
    s_star=0.350,
    zr=1000.0,
    et_max=4.42,
    e_w=0.20,
    interception=2.0,
    name="Prosopis glandulosa",
)

# Table 2. 0.0155 cells/km^2 = 1.55e-8 cells/m^2. a' = 5 km.
RAIN = RainParams(
    storm_rate=0.167,
    mean_cell_depth=25.2,
    cell_density=0.0155 / 1.0e6,
    cell_scale=5000.0,
)

# Initial saturation used in the Fortran input files. Not in the tables.
S_INIT = 0.35

# Equation (14) as printed in Table 4. c1 = 1830 mm/d with c2 = -2.830 mm/d
# produces a loss of hundreds of mm/d between wilting and field capacity,
# which is not physical. c3, c4, c5, and c6 are plausible. Refit before use.
UPSCALED = {
    "c1": 1830.0,
    "c2": -2.830,
    "c3": 1469.0,
    "c4": -781.7,
    "c5": 1.457,
    "c6": 2.689,
}


def mean_storm_depth(rain: RainParams = RAIN) -> float:
    """E[Y] from equation (6), millimeters.

    Table 2 rounds this to 15 mm. With the tabulated λ, a', and E[h] it is
    about 15.3 mm.
    """
    a_km = rain.cell_scale / 1000.0
    density_per_km2 = rain.cell_density * 1.0e6
    return 0.5 * np.pi * density_per_km2 * a_km**2 * rain.mean_cell_depth


def effective_plant(grass_fraction: float, grass: PlantParams = GRASS, woody: PlantParams = WOODY) -> PlantParams:
    """Area-average fluxes and root-zone depth. Volume-average saturation thresholds.

    Matches Table 1. For a 50/50 mix the table gives Zr = 700 mm, ETmax = 4.59 mm/d,
    Ew = 0.165 mm/d, Δ = 1.5 mm, Sw = 0.176, S* = 0.357.
    """
    if not 0.0 <= grass_fraction <= 1.0:
        raise ValueError("grass_fraction must be between 0 and 1")
    fg = float(grass_fraction)
    ft = 1.0 - fg
    zr = fg * grass.zr + ft * woody.zr
    # Storage weights. Porosity is the same for both plants, so it cancels.
    vg = fg * grass.zr / zr
    vt = ft * woody.zr / zr
    return PlantParams(
        s_w=vg * grass.s_w + vt * woody.s_w,
        s_star=vg * grass.s_star + vt * woody.s_star,
        zr=zr,
        et_max=fg * grass.et_max + ft * woody.et_max,
        e_w=fg * grass.e_w + ft * woody.e_w,
        interception=fg * grass.interception + ft * woody.interception,
        name="effective",
    )
