# More sliver fixing

in some instances tiny little grids are tagged along from processing. we
will just give these to an arbitrary nearest feature. They would be
annoying if by chance (... they are about a millionth of the area...) an
arbitrary random point was in them and missed, but the areas themselves
tend to be remarkable inconsequential, and worth randomly reassigning to
a neighbor

## Usage

``` r
reduceFinalGrids(final_grids)
```

## Arguments

- final_grids:

  truthfully nearly final at this point.
