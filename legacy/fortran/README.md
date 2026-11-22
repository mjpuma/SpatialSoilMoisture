# Fortran spatial bucket

This is the spatial bucket model from the last part of Puma (2006). Depths in the input files are centimeters and days. The Python model uses the same physical values in millimeters.

Several working copies of this program existed. This directory is the one kept with the repository. It builds with gfortran.

## Run a short check

From this directory:

```bash
gfortran -std=legacy -O2 -o bucket2d BUCKET2D.F CONTROL.F CUNIF.F EXPO.F \
    INPUT.F LOCATE.F POISSON.F SOILPROP.F STORM.F VEGPROP.F
cd check && ../bucket2d
```

`check/` is a 1 km block for 30 days, so the run finishes quickly. The program still allocates the full sample arrays (about 1 GB). The La Copita case in the parent directory is a 30 km block for 20200 days. That is the dissertation experiment, and it is a long run. Run it from this directory, not from `check/`:

```bash
./bucket2d
```

The program reads `bucket.in` from the current directory. That file names `control.in`, `soil.in`, and `veget.in`.

## What the check should show

The screen prints one line per day: the day, the grass fraction, the next storm day, the mean saturation of the whole block, and the fraction of points above field capacity. Mean saturation stays between 0 and 1.

These files are written in the current directory:

| File | What it is |
|---|---|
| `vegmap.out` | Vegetation type at each sample point. 1 is grass, 2 is woody. Both values appear for the check inputs. |
| `sat.out` | Mean relative saturation on windows that grow from the center of the grid. The last number on each line is the whole block. |
| `et.out`, `leak.out`, `rain.out` | Area-mean evapotranspiration, leakage, and rainfall for those same windows, in cm/d or cm. Evapotranspiration and leakage are the rates at the start of the day, before the soil dries. |
| `varsat.out`, `varet.out`, `varrain.out`, `varleak.out` | Spatial variances on the same windows. |
| `wetf.out` | Fraction of sample points above field capacity. |

A run that stops on the first day, prints saturations outside [0, 1], or leaves `vegmap.out` with only one vegetation type has not exercised the model. The numbers change from run to run because the random sequence is not seeded.

## Input units

`control.in`, one value per line: time step (days), last day, grid width and height (m), plant-grid spacing (m), storms per day, rain cells per m², cell parameter a in 1/m, mean cell-center depth (cm), woody plants per m², mean crown radius (m). The water balance is solved on every 10th plant-grid point. That factor is set in `BUCKET2D.F` as `nspars = 10`.

`soil.in`: initial saturation, hygroscopic saturation, field capacity, porosity, saturated conductivity (cm/d), Clapp-Hornberger b.

`veget.in`: grass, then woody. For each, wilting saturation, stomatal-closure saturation, root depth (cm), maximum evapotranspiration (cm/d), evaporation at wilting (cm/d), interception (cm).

`check/control.in` uses a smaller block and a crown radius large enough that some sample points, which are 100 m apart, fall inside a crown. The parent `control.in` is the La Copita setup: 30 km, mean crown radius 1.5 m, 0.04 woody plants per m².
