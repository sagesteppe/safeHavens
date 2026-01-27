# Identify clusters of populations least separated by landscape resistance

Creates a matrix, and supporting data structures, of a Isolation by
Distance distances using a simplified framework. Given a raster surface
parameterizing resistance to movement (e.g., oceans, lakes, rivers,
topographic roughness, or habitat suitability), this function calculates
pairwise landscape resistance among populations and identifies clusters
of populations that are least isolated from each other.

The function can internally build a resistance raster, but for
performance across multiple species, it is recommended to provide
pre-computed surfaces. Resistance is treated as isotropic (no slope
directionality). Clustering is based on pairwise IBR distances and can
be visualized or used to guide sampling schemes.

## Usage

``` r
populationResistance(
  populations_sf,
  base_raster,
  buffer_dist = 25000,
  planar_proj = NULL,
  n_points = 150,
  resistance_surface = NULL,
  oceans = NULL,
  lakes = NULL,
  rivers = NULL,
  tri = NULL,
  habitat = NULL,
  w_ocean = 100,
  w_lakes = 50,
  w_rivers = 20,
  w_tri = 1,
  w_habitat = 1,
  graph_method = c("complete", "delaunay"),
  ibr_method = "leastcost"
)
```

## Arguments

- populations_sf:

  sf/tibble/data.frame. Coordinates of populations to be clustered. Must
  contain POINT geometries.

- base_raster:

  SpatRaster. Base raster for the study area. Provides template geometry
  and resolution.

- buffer_dist:

  Numeric. Distance to buffer population points to assign likely
  occupied cells (optional; ignored if using `sample_population_cells`
  directly).

- planar_proj:

  Numeric. EPSG code used to project `populations_sf` for buffering
  operations.

- n_points:

  Numeric (defaults to 50). Number of points to sample for creating the
  surface. Increasing points drastically increases compute time.

- resistance_surface:

  SpatRaster. Optional pre-computed resistance raster. If provided, the
  raster-building arguments are ignored.

- oceans:

  SpatRaster. Binary (0/1) raster for ocean cells. Used to increase
  movement cost.

- lakes:

  SpatRaster. Binary (0/1) raster for lakes.

- rivers:

  SpatRaster. Binary (0/1) raster for rivers.

- tri:

  SpatRaster. Continuous raster of topographic roughness (TRI). Used to
  increase cost in mountainous terrain.

- habitat:

  SpatRaster. Continuous raster of habitat suitability. Low values
  increase cost.

- w_ocean:

  Numeric. Weight applied to oceans (default 2000).

- w_lakes:

  Numeric. Weight applied to lakes (default 200).

- w_rivers:

  Numeric. Weight applied to rivers (default 20).

- w_tri:

  Numeric. Weight applied to TRI (default 1).

- w_habitat:

  Numeric. Weight applied to habitat suitability (default 1).

- graph_method:

  Character. One of `"delaunay"` or `"complete"`. Defines the spatial
  graph connecting population points. Delauney produces a sparser graph
  can deal with higher `n`, 'complete' is more computationally intensive
  and n \> 175 take considerable time.

- ibr_method:

  Character. Currently only `"leastcost"`.

## Value

A list with the following elements:

- resistance_raster:

  The parameterized resistance SpatRaster used for IBR calculations.

- sampled_points:

  sf POINT object of population cells used for graph calculations.

- spatial_graph:

  igraph object representing the population connectivity graph.

- edge_list:

  Data.frame of graph edges used for IBR computation.

- ibr_matrix:

  Symmetric numeric matrix of pairwise landscape resistance between
  populations.

## Examples

``` r
if (FALSE) { # \dontrun{
# Prepare resistance raster once

# this also can run internally in `population resistance`, 
# but for time sakes is best to prep ahead of time
# especially if treating multiple species in the same domain. 
res <- buildResistanceSurface(
  base_raster = base_rast,
  oceans = ocean_r,
  lakes = lakes_r,
  rivers = rivers_r,
  tri = tri_r
)

# Run population resistance clustering
out <- populationResistance(
  populations_sf = pops_sf,
  base_raster = base_rast,
  resistance_surface = res,
  n_points = 5,
  graph_method = "mst",
  ibr_method = "leastcost"
)

# Access the results
ibr_matrix <- out$ibr_matrix
graph <- out$spatial_graph
} # }
```
