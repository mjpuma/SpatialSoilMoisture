"""Poisson storms in time and Gaussian rain cells in space. Equations (4)–(7)."""

from __future__ import annotations

import numpy as np

from .parameters import RAIN, RainParams, mean_storm_depth


def spread(distance: np.ndarray, cell_scale: float) -> np.ndarray:
    """g(r) = exp(-2 (r/a')^2). Equation (5)."""
    return np.exp(-2.0 * (np.asarray(distance, dtype=float) / cell_scale) ** 2)


def daily_storm_depths(n_days: int, rng: np.random.Generator, rain: RainParams = RAIN) -> np.ndarray:
    """One exponential storm depth on a day that contains a storm, else zero.

    The probability of a day with at least one storm is 1 - exp(-λ t).
    Depths use E[Y] from equation (6), the spatial mean of the cell model,
    which is the uniform-rainfall case in the paper.
    """
    occurs = rng.random(n_days) < (1.0 - np.exp(-rain.storm_rate))
    depth = rng.exponential(mean_storm_depth(rain), size=n_days)
    return np.where(occurs, depth, 0.0)


def storm_on_point(
    n_storms: int,
    rng: np.random.Generator,
    rain: RainParams = RAIN,
    half_width: float | None = None,
) -> np.ndarray:
    """Rainfall depth at one interior point for ``n_storms`` independent storms.

    Cells are a Poisson process on a square centered on the point. The default
    half-width is 4 a', which retains essentially the whole spread function.
    """
    if half_width is None:
        half_width = 4.0 * rain.cell_scale
    area = (2.0 * half_width) ** 2
    expected = rain.cell_density * area
    out = np.empty(n_storms)
    for i in range(n_storms):
        n_cells = int(rng.poisson(expected))
        if n_cells == 0:
            out[i] = 0.0
            continue
        cx = rng.uniform(-half_width, half_width, n_cells)
        cy = rng.uniform(-half_width, half_width, n_cells)
        depth = rng.exponential(rain.mean_cell_depth, n_cells)
        distance = np.hypot(cx, cy)
        out[i] = np.sum(depth * spread(distance, rain.cell_scale))
    return out


def storm_field(
    x: np.ndarray,
    y: np.ndarray,
    rng: np.random.Generator,
    rain: RainParams = RAIN,
    margin: float | None = None,
) -> np.ndarray:
    """One storm on a rectangular sample grid. ``x`` and ``y`` are meters.

    Cell centers are drawn on the grid extended by ``margin`` so edge points
    are not systematically dry. Returns millimeters with shape (len(y), len(x)).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if margin is None:
        margin = rain.cell_scale
    xmin, xmax = float(x.min()), float(x.max())
    ymin, ymax = float(y.min()), float(y.max())
    area = (xmax - xmin + 2.0 * margin) * (ymax - ymin + 2.0 * margin)
    n_cells = int(rng.poisson(rain.cell_density * area))
    field = np.zeros((y.size, x.size))
    if n_cells == 0:
        return field
    cx = rng.uniform(xmin - margin, xmax + margin, n_cells)
    cy = rng.uniform(ymin - margin, ymax + margin, n_cells)
    depth = rng.exponential(rain.mean_cell_depth, n_cells)
    xx, yy = np.meshgrid(x, y)
    # Chunk so a large cell count does not build one huge distance array.
    for start in range(0, n_cells, 32):
        stop = min(start + 32, n_cells)
        dist = np.hypot(xx[None] - cx[start:stop, None, None], yy[None] - cy[start:stop, None, None])
        field += np.sum(depth[start:stop, None, None] * spread(dist, rain.cell_scale), axis=0)
    return field
