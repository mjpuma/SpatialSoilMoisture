<p align="center"><img src="assets/spatialsoil-logo.png" alt="Spatial Soil logo: rain, grass, a woody plant, and soil" width="620"></p>

# Spatial Soil

Daily soil moisture, evapotranspiration, and leakage when rainfall and vegetation vary in space. The land surface is a mosaic of grass and woody plants. Storms are collections of rain cells. Each plant has its own root zone, and the root zones do not exchange water.

This is the spatial model from Michael J. Puma's Ph.D. dissertation and from the arXiv preprint of that chapter. The Python package follows the equations in that work. The Fortran program that ran the original experiments is in [`legacy/`](legacy/), with the 2004 source kept as it was found and a corrected `BUCKET2D.F` beside it.

## Source

Cite the dissertation. It is the scholarly record, and this model is its last part.

> Puma, M.J. (2006). *Space-time scaling properties of soil moisture and evapotranspiration in water-limited ecosystems.* Ph.D. dissertation, Department of Civil and Environmental Engineering, Princeton University. Advisers: Michael A. Celia and Ignacio Rodriguez-Iturbe.

The dissertation has three parts: rainfall time resolution and the plant-scale water balance, soil moisture averaged over different depths, and this spatial problem. The first two were published as Puma et al. (2006), *Transport in Porous Media*, and Puma et al. (2005), *Advances in Water Resources*. The spatial part was not published in a journal. It was posted in 2016, without the revision that would have followed review, as [arXiv:1606.05256](https://arxiv.org/abs/1606.05256).

Use the dissertation when the citation needs to be the work itself. Use the arXiv when a reader needs a paper they can download. Citing both is the right default: the dissertation for the research, the arXiv for the open text of the spatial chapter. The two are the same study. They are not a substitute for each other if a reader has to find the document.

## What the model does

Depths in the Python code are millimeters and time is days. Parameters are the La Copita sandy loam, *Paspalum setaceum*, *Prosopis glandulosa*, and growing-season thunderstorms used in the dissertation.

- Point water balance: piecewise evapotranspiration, exponential leakage below field capacity, an interception threshold, saturation-excess runoff, and the analytical drying solution of Laio et al. (2001).
- Uniform rain on grass and woody plants together, and the same rain on one column that uses averaged parameters. Saturation is weighted by root-zone storage. Fluxes are area averages.
- Gaussian rain cells and woody crowns placed as a Boolean model of exponential disks.

The 20,000-day runs on a 30 km domain are not repeated here.

With the same rain on both plants, the shallow grass root zone can leak while the storage-weighted mean saturation is still below field capacity. A single effective column, given that mean saturation, does not leak. In a 4,000-day run this happens on 93 days. Leakage is 2.1% of throughfall for the mosaic and 0.6% for the effective column. The daily water balance closes to roundoff. Those figures are the uniform-rain case. The 90% / 8% / 1% split in the dissertation is the case with spatially variable rain, which is not rerun here.

Table 4 of the arXiv preprint cannot be reused as printed. With c1 = 1830 mm/d and c2 = −2.830 mm/d, the loss between wilting and field capacity is about 550 mm/d. In that range leakage is zero and evapotranspiration is a few millimeters per day. A slope of 18.30 mm/d would be plausible. The evapotranspiration coefficients c5 and c6 are plausible as printed.

## Run

```bash
PYTHONPATH=src python -m pytest
PYTHONPATH=src python examples/uniform_rain.py
```

Python 3.10 or later, and numpy. The tests include a comparison of the analytical drying solution with a numerical integration of the same equation.

## Legacy Fortran

[`legacy/fortran/original/`](legacy/fortran/original/) is the August 2004 source as found, including the Silverfrost project file. [`legacy/fortran/BUCKET2D.F`](legacy/fortran/BUCKET2D.F) is that program with the errors corrected: woody crowns are type 2, each storm starts from an empty rainfall field, and spatial means use the row that belongs to the averaging window. Details are in [`legacy/README.md`](legacy/README.md). The input files are in centimeters.

## Relation to tessera

[tessera](https://github.com/mjpuma/tessera) uses the same point-scale loss function, integrated numerically for one root zone shared by several crop varieties. Its rain is a daily marked Poisson process, correlated between farms by a Gaussian copula. This model builds the rain field from explicit cells, which is the construction tessera's documentation points to. Sampling a storm field from this model at farm locations would replace that copula when the question is whether the geometry of rain cells changes the mixture results. That coupling is not built.
