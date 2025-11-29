# place random points in the polygon which will be dissolved with the larger polygons

This function is ran within `GridBasedSample` it will place points
throughout polygon geometries which should be merged to larger polygons
and assign them to neighboring polygons based on how much area we want
to grow these polygons too.

## Usage

``` r
assignGrid_pts(neighb_grid, focal_grid, props, nf_pct)
```

## Arguments

- neighb_grid:
- focal_grid:
- props:
- nf_pct:
