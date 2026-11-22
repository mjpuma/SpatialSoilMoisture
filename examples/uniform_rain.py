"""Uniform rainfall on grass and woody columns, the case in section 4.1.

Run from the repository root:

    PYTHONPATH=src python examples/uniform_rain.py
"""

import numpy as np

from spatial_soil.parameters import S_INIT, SOIL
from spatial_soil.rainfall import daily_storm_depths
from spatial_soil.simulate import simulate_pair


def main() -> None:
    rng = np.random.default_rng(11)
    rain = daily_storm_depths(4000, rng)
    out = simulate_pair(rain, S_INIT)
    leaking = (out["leakage_rate_mean"] > 0.1) & (out["s_mean"] < SOIL.s_fc)
    # Partition of the mosaic, area-weighted.
    rain_total = float(np.sum(rain))
    et = float(np.sum(out["et_mean"]))
    leak = float(np.sum(out["leakage_mean"]))
    runoff = float(0.5 * (np.sum(out["grass"].runoff) + np.sum(out["woody"].runoff)))
    interception = float(0.5 * (np.sum(out["grass"].interception) + np.sum(out["woody"].interception)))
    through = rain_total - interception
    print(f"days: {rain.size}")
    print(f"water-balance residual: {out['balance_error']:.3e} mm")
    print(f"days with leakage while <S> < Sfc: {int(np.count_nonzero(leaking))}")
    print(f"mosaic  ET / leakage / runoff of throughfall: {et / through:.3f} / {leak / through:.3f} / {runoff / through:.3f}")
    eff_through = rain_total - float(np.sum(out["effective"].interception))
    print(
        "effective ET / leakage / runoff of throughfall: "
        f"{out['effective'].et.sum() / eff_through:.3f} / "
        f"{out['effective'].leakage.sum() / eff_through:.3f} / "
        f"{out['effective'].runoff.sum() / eff_through:.3f}"
    )


if __name__ == "__main__":
    main()
