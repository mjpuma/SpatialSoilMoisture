"""Point-scale water balance.

Use `daily_step(s, rain_mm, column)` for one day. The column is a `Column`
of soil and plant parameters, in millimeters and days. The call returns
the saturation at the end of the day, the millimeters of evapotranspiration
and leakage lost while the soil dried, the instantaneous rates at the
start of the day, and the interception and runoff.

The order inside the day is fixed, and it is the order in the dissertation.
Record the rates from the saturation at the start of the day, dry for one
day, then add rain. Interception is a depth threshold. Water above
saturation is runoff.

Check a result with the water balance, in millimeters:

    porosity * Zr * (s_end - s_start) + ET + leakage + interception + runoff - rain

That residual is a roundoff error when the step is right. `tests/test_point.py`
also compares `dry` with a numerical integration of the same loss function.
Passing that test is the check that the analytical drying solution is the
equation that was intended.

Between storms the soil dries by the analytical solution of Laio et al. (2001).
Equations (3) and (8)–(12) of the spatial part of Puma (2006) and of arXiv:1606.05256.
"""

from __future__ import annotations

from typing import NamedTuple, Union

import numpy as np

from .parameters import PlantParams, SoilParams

ArrayLike = Union[float, np.ndarray]


class Column(NamedTuple):
    soil: SoilParams
    plant: PlantParams

    @property
    def capacity(self) -> float:
        """Saturated root-zone storage, mm."""
        return self.soil.porosity * self.plant.zr


class Drying(NamedTuple):
    s: np.ndarray
    et: np.ndarray
    leakage: np.ndarray


class Day(NamedTuple):
    s: np.ndarray
    et: np.ndarray
    leakage: np.ndarray
    et_rate: np.ndarray
    leakage_rate: np.ndarray
    interception: np.ndarray
    runoff: np.ndarray


def _as_array(s: ArrayLike) -> tuple[np.ndarray, bool]:
    arr = np.asarray(s, dtype=float)
    return arr, arr.ndim == 0


def et_rate(s: ArrayLike, column: Column) -> np.ndarray:
    """Instantaneous evapotranspiration, mm/d. Equation (10)."""
    p = column.plant
    s_h = column.soil.s_h
    s_arr, scalar = _as_array(s)
    et = np.zeros_like(s_arr, dtype=float)
    evaporating = (s_arr > s_h) & (s_arr <= p.s_w)
    stressed = (s_arr > p.s_w) & (s_arr <= p.s_star)
    unlimited = s_arr > p.s_star
    et[evaporating] = p.e_w * (s_arr[evaporating] - s_h) / (p.s_w - s_h)
    et[stressed] = p.e_w + (p.et_max - p.e_w) * (s_arr[stressed] - p.s_w) / (p.s_star - p.s_w)
    et[unlimited] = p.et_max
    if scalar:
        return et.reshape(())
    return et


def leakage_rate(s: ArrayLike, column: Column) -> np.ndarray:
    """Instantaneous leakage, mm/d. Equation (11). Zero at and below field capacity."""
    soil = column.soil
    s_arr, scalar = _as_array(s)
    leak = np.zeros_like(s_arr, dtype=float)
    wet = s_arr > soil.s_fc
    if np.any(wet):
        beta = 2.0 * soil.b + 4.0
        scale = np.exp(beta * (1.0 - soil.s_fc)) - 1.0
        leak[wet] = soil.ksat / scale * (np.exp(beta * (s_arr[wet] - soil.s_fc)) - 1.0)
    if scalar:
        return leak.reshape(())
    return leak


def _check_order(column: Column) -> None:
    soil, plant = column.soil, column.plant
    if not soil.s_h < plant.s_w < plant.s_star < soil.s_fc < 1.0:
        raise ValueError("need s_h < s_w < s_star < s_fc < 1 for the piecewise solution")


def dry(s: ArrayLike, dt: float, column: Column) -> Drying:
    """Analytical drying over ``dt`` days. No rainfall in this interval."""
    if dt < 0.0:
        raise ValueError("dt must be non-negative")
    _check_order(column)
    s_arr, scalar = _as_array(s)
    shape = s_arr.shape
    s_now = np.array(s_arr, dtype=float, copy=True).reshape(-1)
    et = np.zeros_like(s_now)
    leak = np.zeros_like(s_now)
    remaining = np.full(s_now.shape, float(dt))

    soil, plant = column.soil, column.plant
    cap = column.capacity
    eta = plant.et_max / cap
    eta_w = plant.e_w / cap
    beta = 2.0 * soil.b + 4.0
    m = (soil.ksat / cap) / (np.exp(beta * (1.0 - soil.s_fc)) - 1.0)
    a = eta - m  # coefficient in ds/dt = -a - m exp(beta (s-s_fc))

    # Above field capacity: ET at ETmax and leakage.
    idx = np.flatnonzero((s_now > soil.s_fc) & (remaining > 0.0))
    if idx.size and abs(a) > 0.0:
        ss = s_now[idx]
        temp = m * np.exp(beta * (ss - soil.s_fc))
        t_fc = (beta * (ss - soil.s_fc) + np.log(eta / (a + temp))) / (a * beta)
        t_fc = np.maximum(t_fc, 0.0)
        t_spend = np.minimum(remaining[idx], t_fc)
        numer = (a + temp) * np.exp(a * beta * t_spend) - temp
        s_new = ss - np.log(numer / a) / beta
        storage_loss = cap * (ss - s_new)
        et_piece = plant.et_max * t_spend
        et[idx] += et_piece
        leak[idx] += storage_loss - et_piece
        crossed = remaining[idx] >= t_fc
        s_now[idx] = np.where(crossed, soil.s_fc, s_new)
        remaining[idx] = np.where(crossed, remaining[idx] - t_fc, 0.0)

    # Field capacity down to incipient stomatal closure: ET = ETmax, no leakage.
    idx = np.flatnonzero((s_now > plant.s_star) & (remaining > 0.0))
    if idx.size:
        t_star = (s_now[idx] - plant.s_star) / eta
        t_spend = np.minimum(remaining[idx], t_star)
        s_new = s_now[idx] - eta * t_spend
        et[idx] += cap * (s_now[idx] - s_new)
        crossed = remaining[idx] >= t_star
        s_now[idx] = np.where(crossed, plant.s_star, s_new)
        remaining[idx] = np.where(crossed, remaining[idx] - t_star, 0.0)

    # Stressed transpiration, down to wilting.
    bb = (eta - eta_w) / (plant.s_star - plant.s_w)
    idx = np.flatnonzero((s_now > plant.s_w) & (remaining > 0.0))
    if idx.size:
        rate0 = eta_w + bb * (s_now[idx] - plant.s_w)
        t_w = np.log(rate0 / eta_w) / bb
        t_spend = np.minimum(remaining[idx], t_w)
        s_new = (plant.s_w - eta_w / bb) * (1.0 - np.exp(-bb * t_spend)) + s_now[idx] * np.exp(-bb * t_spend)
        et[idx] += cap * (s_now[idx] - s_new)
        crossed = remaining[idx] >= t_w
        s_now[idx] = np.where(crossed, plant.s_w, s_new)
        remaining[idx] = np.where(crossed, remaining[idx] - t_w, 0.0)

    # Evaporation below wilting, approaching the hygroscopic point.
    gamma = eta_w / (plant.s_w - soil.s_h)
    idx = np.flatnonzero((s_now > soil.s_h) & (remaining > 0.0))
    if idx.size:
        s_new = soil.s_h + (s_now[idx] - soil.s_h) * np.exp(-gamma * remaining[idx])
        et[idx] += cap * (s_now[idx] - s_new)
        s_now[idx] = s_new
        remaining[idx] = 0.0

    if scalar:
        return Drying(s_now.reshape(()), et.reshape(()), leak.reshape(()))
    return Drying(s_now.reshape(shape), et.reshape(shape), leak.reshape(shape))


def apply_rain(s: ArrayLike, rain: ArrayLike, column: Column) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interception threshold, then saturation excess. Equations (8) and (9).

    Returns the new saturation, intercepted depth, and runoff depth.
    """
    s_arr = np.asarray(s, dtype=float)
    rain_arr = np.asarray(rain, dtype=float)
    intercepted = np.minimum(rain_arr, column.plant.interception)
    through = rain_arr - intercepted
    trial = s_arr + through / column.capacity
    runoff = np.maximum(trial - 1.0, 0.0) * column.capacity
    s_new = np.minimum(trial, 1.0)
    return s_new, intercepted, runoff


def daily_step(s: ArrayLike, rain: ArrayLike, column: Column) -> Day:
    """One day: record instantaneous fluxes, dry for one day, then add rain."""
    rates_et = et_rate(s, column)
    rates_l = leakage_rate(s, column)
    dried = dry(s, 1.0, column)
    s_new, intercepted, runoff = apply_rain(dried.s, rain, column)
    return Day(
        s=s_new,
        et=dried.et,
        leakage=dried.leakage,
        et_rate=rates_et,
        leakage_rate=rates_l,
        interception=intercepted,
        runoff=runoff,
    )
