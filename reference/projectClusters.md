# Project current environmental clusters onto a future climate scenario

Main entry point for the future-projection workflow. The steps are:

## Usage

``` r
projectClusters(
  eSDM_object,
  current_clusters,
  future_predictors,
  current_predictors,
  planar_proj,
  coord_wt = 0.001,
  mess_threshold = 0,
  cluster_novel = TRUE,
  n_novel_pts = 500,
  n_sample_per_cluster = 50,
  nbclust_args = list(),
  thresh_metric = "sensitivity",
  thresholds
)
```

## Arguments

- eSDM_object:

  output from
  [`elasticSDM()`](https://sagesteppe.github.io/safeHavens/reference/elasticSDM.md).

- current_clusters:

  Output list from
  [`EnvironmentalBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EnvironmentalBasedSample.md).

- future_predictors:

  SpatRaster of future climate. Layer names must match those retained in
  `current_model`.

- current_predictors:

  SpatRaster of current climate (used for MESS reference and for
  `rescaleFuture` standardisation).

- planar_proj:

  EPSG code or proj4 string for a planar projection in metres (same as
  used in the current analysis).

- coord_wt:

  Numeric, default `0.001`. Coordinate weighting passed to
  `add_weighted_coordinates`.

- mess_threshold:

  MESS values below this are treated as novel climate. Default `0`.

- cluster_novel:

  Boolean, defualt `TRUE`. If `TRUE` and novel cells exist, cluster them
  independently with NbClust. If `FALSE` novel cells are left as `NA`.

- n_novel_pts:

  Numeric, default 500. Number of points to sample from novel areas for
  clustering.

- n_sample_per_cluster:

  Number of points to sample from each cluster (existing *and* novel)
  for the relationship tree. Default `50`.

- nbclust_args:

  Named list of arguments forwarded to
  [`NbClust::NbClust`](https://rdrr.io/pkg/NbClust/man/NbClust.html).
  Sensible defaults are set internally (`min.nc = 2`, `max.nc = 10`,
  `method = "ward.D2"`, `index = "all"`).

- thresh_metric:

  Character. Default 'sensitivity', dismo::threshold to use for cutting
  the future sdm into a binary surface.

- thresholds:

  current era thresholds from `postProcessSDM`

## Value

A named list:

- clusters_sf:

  `sf` polygons with column `ID`.

- suitable_habitat:

  Raster of masked suitable habitat under future conditions

- novel_mask:

  `SpatRaster` — logical, `TRUE` where MESS \< threshold.

- mess:

  `SpatRaster` — raw MESS scores (minimum across all variables).

- changes:

  `data.frame` — per-cluster area and centroid-shift metrics.

- novel_similarity:

  `data.frame` — nearest existing cluster and mean silhouette width for
  each novel cluster. Zero-row data frame if none.

## Details

1.  Rescale future predictors with current betas (`rescaleFuture`).

2.  Run [`dismo::mess`](https://rdrr.io/pkg/dismo/man/mess.html) to
    identify novel climate cells.

3.  Predict known-climate cells with the existing KNN classifier.

4.  If novel cells exist and `cluster_novel = TRUE`, cluster them
    independently with `NbClust` (`cluster_novel_areas`).

5.  Sample points from every cluster (existing + novel), force them onto
    a single tree, and extract nearest-existing-cluster relationships
    via silhouette (`analyze_cluster_relationships`).

6.  Polygonise and calculate area / centroid changes.
