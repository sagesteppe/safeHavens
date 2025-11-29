# turn the point grid suggestions made by assignGrid_pts into polygons

This function is part of the grid based sampling process to turn small
grid cells, which are to be broken up, into a larger existing grid
cells.

## Usage

``` r
snapGrids(x, neighb_grid, focal_grid)
```

## Arguments

- x:

  output of `assignGrid_pts`

- neighb_grid:

  the neighboring grid options.

- focal_grid:

  the grid to reassign the area of.
