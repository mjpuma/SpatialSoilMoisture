"""Checks for one root zone.

`test_analytical_drying_matches_rk4` integrates the loss function
independently and compares it with `dry`. A failure there means the
closed-form drying solution has drifted from the equation.

`test_daily_water_balance_with_rain_and_runoff` checks

    storage change + ET + leakage + interception + runoff - rain = 0

on a day that includes runoff. The residual must be a roundoff error.
"""

import numpy as np

from spatial_soil.parameters import GRASS, SOIL, WOODY
from spatial_soil.point import Column, apply_rain, daily_step, dry, et_rate, leakage_rate


def rk4(s, dt, column, n=4000):
    h = dt / n
    cap = column.capacity

    def f(value):
        return -float(et_rate(value, column) + leakage_rate(value, column)) / cap

    for _ in range(n):
        k1 = f(s)
        k2 = f(s + 0.5 * h * k1)
        k3 = f(s + 0.5 * h * k2)
        k4 = f(s + h * k3)
        s = s + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return s


def test_et_breakpoints():
    column = Column(SOIL, GRASS)
    assert et_rate(0.10, column) == 0.0
    assert et_rate(SOIL.s_h, column) == 0.0
    assert np.isclose(et_rate(GRASS.s_w, column), GRASS.e_w)
    assert np.isclose(et_rate(GRASS.s_star, column), GRASS.et_max)
    assert et_rate(0.9, column) == GRASS.et_max


def test_leakage_endpoints():
    column = Column(SOIL, GRASS)
    assert leakage_rate(SOIL.s_fc, column) == 0.0
    assert leakage_rate(0.3, column) == 0.0
    assert np.isclose(leakage_rate(1.0, column), SOIL.ksat)


def test_analytical_drying_matches_rk4():
    for plant in (GRASS, WOODY):
        column = Column(SOIL, plant)
        for s0 in (0.95, 0.75, 0.50, 0.30, 0.20, 0.15):
            for dt in (0.2, 1.0, 10.0):
                analytical = float(dry(s0, dt, column).s)
                numerical = rk4(s0, dt, column)
                assert abs(analytical - numerical) < 2e-3, (plant.name, s0, dt, analytical, numerical)


def test_drying_conserves_storage_as_et_plus_leakage():
    column = Column(SOIL, GRASS)
    s0 = 0.9
    dried = dry(s0, 1.0, column)
    storage_loss = column.capacity * (s0 - float(dried.s))
    assert abs(storage_loss - float(dried.et + dried.leakage)) < 1e-8


def test_no_loss_at_hygroscopic_point():
    column = Column(SOIL, WOODY)
    dried = dry(SOIL.s_h, 5.0, column)
    assert float(dried.s) == SOIL.s_h
    assert float(dried.et) == 0.0
    assert float(dried.leakage) == 0.0


def test_daily_water_balance_with_rain_and_runoff():
    column = Column(SOIL, GRASS)
    s0 = 0.8
    rain = 80.0  # enough to hit saturation after drying
    day = daily_step(s0, rain, column)
    # Reconstruct the pre-rain state from the integrated losses.
    s_after_dry = s0 - float(day.et + day.leakage) / column.capacity
    s_new, intercepted, runoff = apply_rain(s_after_dry, rain, column)
    assert np.isclose(day.s, s_new)
    residual = column.capacity * (float(day.s) - s0) + float(day.et + day.leakage + day.interception + day.runoff) - rain
    assert abs(residual) < 1e-8
    assert float(day.runoff) > 0.0
    assert float(day.s) <= 1.0


def test_small_rain_is_fully_intercepted():
    column = Column(SOIL, GRASS)
    s_new, intercepted, runoff = apply_rain(0.4, 0.4, column)
    assert intercepted == 0.4
    assert runoff == 0.0
    assert s_new == 0.4
