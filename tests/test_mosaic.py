"""The uniform-rain, two-vegetation result in section 4.1."""

import numpy as np

from spatial_soil.parameters import GRASS, SOIL, S_INIT, WOODY, effective_plant
from spatial_soil.point import Column, leakage_rate
from spatial_soil.rainfall import daily_storm_depths
from spatial_soil.simulate import simulate_pair


def test_effective_parameters_match_table_1():
    plant = effective_plant(0.5)
    assert plant.zr == 700.0
    assert np.isclose(plant.et_max, 4.59)
    assert np.isclose(plant.e_w, 0.165)
    assert np.isclose(plant.interception, 1.5)
    assert np.isclose(plant.s_w, 0.176, atol=0.001)
    # Table prints 0.357. Volume weighting gives 0.356.
    assert np.isclose(plant.s_star, 0.357, atol=0.002)

    # 56.6% grass is the Boolean cover implied by λ = 0.04 m^-2 and μ = 1.5 m.
    savanna = effective_plant(0.566)
    assert np.isclose(savanna.zr, 660.4, atol=0.2)
    assert np.isclose(savanna.et_max, 4.61, atol=0.01)


def test_leakage_occurs_below_mean_field_capacity():
    """A grass plant can leak while the storage-weighted mean saturation is below Sfc.

    Grass at 0.70 and woody at 0.40:
    <S> = (400*0.70 + 1000*0.40) / 1400 = 0.486, below Sfc = 0.56.
    Only the grass is above field capacity, so area-averaged leakage is positive.
    The effective column, given that same <S>, does not leak.
    """
    s_grass = 0.70
    s_woody = 0.40
    s_mean = (GRASS.zr * s_grass + WOODY.zr * s_woody) / (GRASS.zr + WOODY.zr)
    assert s_mean < SOIL.s_fc
    leak = 0.5 * (
        float(leakage_rate(s_grass, Column(SOIL, GRASS))) + float(leakage_rate(s_woody, Column(SOIL, WOODY)))
    )
    assert leak > 0.0
    assert float(leakage_rate(s_mean, Column(SOIL, effective_plant(0.5)))) == 0.0


def test_pair_simulation_closes_and_leaks_below_field_capacity():
    rng = np.random.default_rng(11)
    rain = daily_storm_depths(4000, rng)
    out = simulate_pair(rain, S_INIT)
    assert out["balance_error"] < 1e-6
    leaking_while_dry = (out["leakage_rate_mean"] > 0.1) & (out["s_mean"] < SOIL.s_fc)
    assert np.count_nonzero(leaking_while_dry) > 20
    # Over a long run the mosaic leaks more than the single effective column.
    assert out["leakage_mean"].sum() > out["effective"].leakage.sum()
