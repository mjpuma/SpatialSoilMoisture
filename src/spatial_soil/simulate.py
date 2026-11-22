"""Daily simulations of one column, of a grass–woody pair, and of the effective column.

`simulate_column(rain, column, s_init)` integrates one root zone.
`rain` is a 1-D array of daily depths in millimeters. The function returns
the series and the daily water-balance residual. A residual larger than
about 1e-6 mm means the day did not conserve water.

`simulate_pair(rain, s_init)` puts the same rain on grass and on the woody
plant, and on a third column that uses the averaged parameters. The
returned `s_mean` weights the two saturations by root-zone storage.
`et_mean` and `leakage_mean` are area means of the millimeter totals.
`leakage_rate_mean` is the area mean of the instantaneous leakage at the
start of each day.

The dissertation result to look for: some days have `leakage_rate_mean > 0`
while `s_mean` is still below field capacity, and the mosaic leaks more
over the run than the effective column. `examples/uniform_rain.py` prints
that comparison. `tests/test_mosaic.py` asserts it.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from .averaging import area_mean, saturation_mean
from .parameters import GRASS, SOIL, WOODY, PlantParams, SoilParams, effective_plant
from .point import Column, apply_rain, dry, et_rate, leakage_rate


class Series(NamedTuple):
    s: np.ndarray
    et: np.ndarray
    leakage: np.ndarray
    et_rate: np.ndarray
    leakage_rate: np.ndarray
    runoff: np.ndarray
    interception: np.ndarray
    rain: np.ndarray


def simulate_column(rain: np.ndarray, column: Column, s_init: float) -> tuple[Series, np.ndarray]:
    """Integrate one column. Returns the series and the daily water-balance residual, mm.

    Residual is φ Zr (s_after - s_before) + ET + L + interception + runoff - rain.
    """
    rain = np.asarray(rain, dtype=float)
    n = rain.shape[0]
    s = np.empty(n)
    et = np.empty(n)
    leak = np.empty(n)
    et_r = np.empty(n)
    leak_r = np.empty(n)
    runoff = np.empty(n)
    interception = np.empty(n)
    residual = np.empty(n)
    state = float(s_init)
    for t in range(n):
        s_before = state
        et_r[t] = float(et_rate(state, column))
        leak_r[t] = float(leakage_rate(state, column))
        dried = dry(state, 1.0, column)
        state_f, inter_f, run_f = apply_rain(dried.s, rain[t], column)
        state = float(state_f)
        et[t] = float(dried.et)
        leak[t] = float(dried.leakage)
        interception[t] = float(inter_f)
        runoff[t] = float(run_f)
        s[t] = state
        residual[t] = (
            column.capacity * (state - s_before) + et[t] + leak[t] + interception[t] + runoff[t] - rain[t]
        )
    series = Series(s, et, leak, et_r, leak_r, runoff, interception, rain)
    return series, residual


def simulate_pair(
    rain: np.ndarray,
    s_init: float,
    soil: SoilParams = SOIL,
    grass: PlantParams = GRASS,
    woody: PlantParams = WOODY,
    grass_fraction: float = 0.5,
) -> dict:
    """Equal rainfall on grass and woody columns, plus the effective-parameter column.

    This is section 4.1: uniform rain, two vegetation types. Spatial means use
    area weights for fluxes and storage weights for saturation.
    """
    grass_col = Column(soil, grass)
    woody_col = Column(soil, woody)
    eff_col = Column(soil, effective_plant(grass_fraction, grass, woody))
    g, g_err = simulate_column(rain, grass_col, s_init)
    w, w_err = simulate_column(rain, woody_col, s_init)
    e, e_err = simulate_column(rain, eff_col, s_init)

    fg = grass_fraction
    ft = 1.0 - fg
    # Storage per unit area is φ Zr, and φ is common.
    wg = fg * grass.zr
    wt = ft * woody.zr
    s_mean = (wg * g.s + wt * w.s) / (wg + wt)
    et_mean = fg * g.et + ft * w.et
    leak_mean = fg * g.leakage + ft * w.leakage
    et_rate_mean = fg * g.et_rate + ft * w.et_rate
    leak_rate_mean = fg * g.leakage_rate + ft * w.leakage_rate
    return {
        "grass": g,
        "woody": w,
        "effective": e,
        "s_mean": s_mean,
        "et_mean": et_mean,
        "leakage_mean": leak_mean,
        "et_rate_mean": et_rate_mean,
        "leakage_rate_mean": leak_rate_mean,
        "balance_error": np.maximum.reduce([np.max(np.abs(g_err)), np.max(np.abs(w_err)), np.max(np.abs(e_err))]),
        "effective_column": eff_col,
    }


def grid_means(s: np.ndarray, woody: np.ndarray, grass: PlantParams = GRASS, woody_plant: PlantParams = WOODY) -> float:
    zr = np.where(woody, woody_plant.zr, grass.zr)
    return saturation_mean(s, zr)


def flux_mean(field: np.ndarray) -> float:
    return area_mean(field)
