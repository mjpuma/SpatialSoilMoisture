"""Woody crowns on a grass matrix.

`woody_mask(x, y, rng)` returns a boolean array, True where a sample point
lies inside a crown. `x` and `y` are 1-D coordinates in meters. Crowns are
disks whose centers are a Poisson process and whose radii are exponential.
The default intensity, 0.04 per m², and mean radius, 1.5 m, are the
savanna in the dissertation. The cover fraction of that savanna is

    woody_cover_probability(0.04, 1.5)

which is 0.432. A large sample from `woody_mask` should come out near that
value. `tests/test_rainfall_vegetation.py` is that check. The dissertation
text reports 43.4 percent woody cover from a realization of the same model.
"""

from __future__ import annotations

import numpy as np


def woody_cover_probability(density: float, mean_radius: float) -> float:
    """Cover fraction 1 - exp(-λ E[π R^2]). For an exponential radius, E[R^2] = 2 μ^2."""
    mean_area = np.pi * 2.0 * mean_radius**2
    return float(1.0 - np.exp(-density * mean_area))


def woody_mask(
    x: np.ndarray,
    y: np.ndarray,
    rng: np.random.Generator,
    density: float = 0.04,
    mean_radius: float = 1.5,
    margin_factor: float = 12.0,
) -> np.ndarray:
    """Boolean woody cover on the grid defined by 1-D coordinates ``x`` and ``y``.

    Trees are generated in a margin around the grid. The margin is
    ``margin_factor`` mean radii. The exponential tail beyond that is negligible.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    margin = margin_factor * mean_radius
    xmin, xmax = float(x.min()) - margin, float(x.max()) + margin
    ymin, ymax = float(y.min()) - margin, float(y.max()) + margin
    n_trees = int(rng.poisson(density * (xmax - xmin) * (ymax - ymin)))
    cover = np.zeros((y.size, x.size), dtype=bool)
    if n_trees == 0:
        return cover
    cx = rng.uniform(xmin, xmax, n_trees)
    cy = rng.uniform(ymin, ymax, n_trees)
    radius = rng.exponential(mean_radius, n_trees)
    xx, yy = np.meshgrid(x, y)
    for start in range(0, n_trees, 64):
        stop = min(start + 64, n_trees)
        dist2 = (xx[None] - cx[start:stop, None, None]) ** 2 + (yy[None] - cy[start:stop, None, None]) ** 2
        cover |= np.any(dist2 <= radius[start:stop, None, None] ** 2, axis=0)
    return cover
