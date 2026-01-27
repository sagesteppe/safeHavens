# Sample a species based on Isolation by Resistance Distance (IBR)

Create `n` seed collection areas based on the habitat, and geographic,
resistance between points.

## Usage

``` r
IBRSurface(
  base_raster,
  pop_raster,
  resistance_surface,
  pts_sf,
  ibr_matrix,
  fixedClusters = FALSE,
  n = NULL,
  min.nc = 5,
  max.nc = 20,
  planar_proj,
  max_geo_cells = NULL,
  geo_iters = 10,
  boundary_width = 2,
  distance_method = c("haversine", "cosine")
)
```

## Arguments

- base_raster:

  a Raster surface to sample the points within, ....

- pop_raster:

  Spatiall buffered population raster from `populationResistance`.

- resistance_surface:

  a resistance surface raster from `build_resistance_surface` or
  `populationResistance`.

- pts_sf:

  sf/tibble/dataframe of point locations from `populationResistance`
  \$pts_sf

- ibr_matrix:

  Distance matrix, from `populationResistance` slot \$Rmat

- fixedClusters:

  Boolean. Defaults to FALSE, which will create n clusters. If False
  then use NbClust::NbClust to determine the optimal number of clusters.
  TRUE not yet suported.

- n:

  Numeric. the number of clusters desired.

- min.nc:

  Numeric. Minimum number of clusters to test if fixedClusters=FALSE,
  defaults to 5.

- max.nc:

  Numeric. Maximum number of clusters to test if fixedClusters=FALSE,
  defaults to 20.

- planar_proj:

  For numbering clusters in return object.

- max_geo_cells:

  Numeric. Maximum number of cells to grow cheap geographic front growth
  of cluster centers by.

- geo_iters:

  Numeric. Maximum number of iterations to grow hulls by.

- boundary_width:

  Numeric. Maximum boundary width before switching to final assignment
  method.

- distance_method:

  Great circle distance calculation method passed onto ot terra. One of
  'haversine' (default), or 'cosine'.

## Examples

``` r
 # see package vignette
```
