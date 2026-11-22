# Legacy Fortran

`fortran/` is the spatial bucket model written by M. J. Puma in August 2004 (Silverfrost FTN95, project file `original/sparbucket.ftn95p`). It is the program behind the last part of Puma (2006) and arXiv:1606.05256. The copy found in 2026 is the only one located; an earlier or later source may have existed.

`fortran/original/` is that copy, unchanged. `fortran/BUCKET2D.F` is the same program with five corrections. The other source files are unchanged. Compiled objects and the Windows executable are not included.

The input files are in centimeters and days. `Zr = 40` is 40 cm, which is the 400 mm grass root depth in the dissertation and the arXiv paper. `Ksat = 82.2` is 82.2 cm/d, which is 822 mm/d. The Python model uses millimeters.

## Corrections in `BUCKET2D.F`

1. A cell inside a woody crown is vegetation type 2. The 2004 line was `iveggrid(i,k) = 1!2`. In Fortran `!` starts a comment, so the `2` was never compiled and every cell stayed grass.
2. Each storm zeroes the rainfall field before adding rain cells. The 2004 code added a storm on top of the previous one, so two storms on consecutive days were counted twice. The storm generator can produce a one-day gap.
3. Nested spatial means use the first row of the averaging window. The 2004 code read row 1 of the whole grid for that row. The full-domain window starts on row 1, so it was unaffected. Every smaller window was not.
4. A tree that falls on the south edge uses the same crown window as a tree on the east edge.
5. The wet fraction printed to the screen is the fraction for the current day. The file `wetf.out` already was.

`random` is the Silverfrost FTN95 intrinsic. gfortran does not provide it, and the generator is not seeded, so two runs do not match. The daily time step is still one shared variable, reset to 1 after each cell.

Runoff is applied by holding saturation at 1 and is not written out. Evapotranspiration and leakage written by the program are the instantaneous rates at the start of the day, which is what section 2.2.6 of the paper states. The Python model also returns the millimeter totals for the day, which the water balance needs.
