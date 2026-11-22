"""Rainfall moments and the Boolean woody cover."""

import numpy as np

from spatial_soil.averaging import nested_means, saturation_mean
from spatial_soil.parameters import RAIN, mean_storm_depth
from spatial_soil.rainfall import storm_on_point
from spatial_soil.vegetation import woody_cover_probability, woody_mask


def test_table_mean_depth_is_about_15_mm():
    # Equation (6) with Table 2. The table rounds 15.3 mm to 15 mm.
    assert abs(mean_storm_depth(RAIN) - 15.34) < 0.05


def test_point_rainfall_matches_moment_equations():
    rng = np.random.default_rng(7)
    n = 6000
    samples = storm_on_point(n, rng, RAIN)
    mu = mean_storm_depth(RAIN)
    # Equation (7): Var = (1/2) π λ a^2 (E[h])^2 = E[Y] * E[h] for this model.
    variance = mu * RAIN.mean_cell_depth
    assert abs(samples.mean() - mu) < 1.0
    assert abs(samples.var(ddof=1) - variance) / variance < 0.08


def test_woody_cover_matches_boolean_model():
    density = 0.04
    radius = 1.5
    expected = woody_cover_probability(density, radius)
    # The paper reports 43.4% woody cover. The Boolean mean is 43.2%.
    assert abs(expected - 0.432) < 0.002
    rng = np.random.default_rng(3)
    # Spacing is larger than a crown, so sample points are nearly independent.
    x = np.linspace(0.0, 400.0, 80)
    y = np.linspace(0.0, 400.0, 80)
    mask = woody_mask(x, y, rng, density=density, mean_radius=radius)
    assert abs(mask.mean() - expected) < 0.03


def test_saturation_mean_weights_by_root_depth():
    s = np.array([1.0, 0.0])
    zr = np.array([400.0, 1000.0])
    assert saturation_mean(s, zr) == 400.0 / 1400.0


def test_nested_means_follow_the_centered_windows():
    field = np.arange(16, dtype=float).reshape(4, 4)
    means = nested_means(field)
    # Even grid: central 2x2, then the full 4x4.
    center = field[1:3, 1:3]
    assert np.isclose(means[0], center.mean())
    assert np.isclose(means[1], field.mean())
    assert means.shape == (2,)
