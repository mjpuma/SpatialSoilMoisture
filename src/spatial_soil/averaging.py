"""Spatial averages that keep the physical meaning of each variable.

Evapotranspiration, leakage, and rainfall are depths per area, so they are
area averages, equation (1). Relative saturation is water per void volume, so
it is averaged with weights φ Zr A, equation (2). On a uniform grid with one
soil, the weights reduce to root-zone depth.
"""

from __future__ import annotations

import numpy as np


def area_mean(field: np.ndarray, weights: np.ndarray | None = None) -> float:
    field = np.asarray(field, dtype=float)
    if weights is None:
        return float(np.mean(field))
    w = np.asarray(weights, dtype=float)
    return float(np.sum(w * field) / np.sum(w))


def saturation_mean(s: np.ndarray, zr: np.ndarray, porosity: float = 1.0) -> float:
    """Volume-weighted mean relative saturation. Equation (2)."""
    s = np.asarray(s, dtype=float)
    weight = porosity * np.asarray(zr, dtype=float)
    return float(np.sum(weight * s) / np.sum(weight))


def nested_windows(n: int):
    """Square windows centered on an n by n grid, from the center out to the full grid.

    Even n starts from the central 2 by 2 block, matching the averaging areas
    used to study scale in the paper. Yields (lo, hi) slices.
    """
    if n < 1:
        raise ValueError("n must be positive")
    if n % 2 == 0:
        lo = n // 2 - 1
        hi = n // 2 + 1
    else:
        lo = n // 2
        hi = lo + 1
    while lo >= 0 and hi <= n:
        yield lo, hi
        lo -= 1
        hi += 1


def nested_means(field: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    """Mean of ``field`` on each window from ``nested_windows``."""
    field = np.asarray(field, dtype=float)
    if field.ndim != 2 or field.shape[0] != field.shape[1]:
        raise ValueError("nested means are defined on a square field")
    if weights is None:
        weights = np.ones_like(field)
    else:
        weights = np.asarray(weights, dtype=float)
    means = []
    for lo, hi in nested_windows(field.shape[0]):
        w = weights[lo:hi, lo:hi]
        means.append(np.sum(w * field[lo:hi, lo:hi]) / np.sum(w))
    return np.asarray(means)
