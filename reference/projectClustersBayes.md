# Project posterior-clustered environmental space onto a future climate scenario

Bayesian analogue of
[`projectClusters()`](https://sagesteppe.github.io/safeHavens/reference/projectClusters.md).
The steps are:

## Usage

``` r
projectClustersBayes(
  bSDM_object,
  posterior_clusters,
  future_predictors,
  current_predictors,
  threshold_rasts,
  planar_proj,
  coord_wt = NULL,
  mess_threshold = 0,
  cluster_novel = TRUE,
  n_novel_pts = 500,
  n_sample_per_cluster = 50,
  n_future_draws = NULL,
  n_future_pts = 500,
  nbclust_args = list(),
  thresh_metric = "sensitivity"
)
```

## Arguments

- bSDM_object:

  Output list from
  [`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md).

- posterior_clusters:

  Output list from
  [`PosteriorCluster()`](https://sagesteppe.github.io/safeHavens/reference/PosteriorCluster.md).

- future_predictors:

  `SpatRaster` of future climate. Layer names must match those used in
  `bSDM_object`.

- current_predictors:

  `SpatRaster` of current climate (used for MESS reference and MESS
  computation (raw units required).

- threshold_rasts:

  Output list from
  [`PostProcessSDM()`](https://sagesteppe.github.io/safeHavens/reference/PostProcessSDM.md),
  containing: `$FinalRasters` (a `SpatRaster` stack including a binary
  threshold layer) and `$Threshold` (a one-row `data.frame` of threshold
  metric values).

- planar_proj:

  EPSG code or proj4 string for a planar projection in metres (same as
  used in the current analysis).

- coord_wt:

  Ignored. Coordinate weighting is read directly from
  `posterior_clusters$ScalingParams$coord_wt` to guarantee the future
  feature space matches the scale the KNNs were trained on. The
  parameter is retained in the signature for backwards compatibility
  only.

- mess_threshold:

  MESS values below this are treated as novel climate. Default `0`.

- cluster_novel:

  Logical, default `TRUE`. If `TRUE` and novel cells exist, cluster them
  independently with NbClust + KNN. If `FALSE` novel cells are left as
  `NA`.

- n_novel_pts:

  Numeric, default `500`. Number of points to sample from novel areas
  for independent clustering.

- n_sample_per_cluster:

  Number of points to sample from each cluster (existing *and* novel)
  for the silhouette relationship analysis. Default `50`.

- n_future_draws:

  Integer. Number of posterior beta draws to use when computing future
  cluster-assignment stability. Defaults to the number of draws stored
  in `posterior_clusters$ScalingParams$beta_draws` (i.e. however many
  were used in
  [`PosteriorCluster()`](https://sagesteppe.github.io/safeHavens/reference/PosteriorCluster.md)),
  capped at `100`. Pass an explicit integer to override.

- n_future_pts:

  Integer. Number of points to sample from the future suitable-habitat
  surface for the draw loop. These are analogous to the fixed sample
  points in
  [`PosteriorCluster()`](https://sagesteppe.github.io/safeHavens/reference/PosteriorCluster.md)
  and serve the same purpose: stability is computed on this point set,
  then painted to raster once via KNN. Default `500`.

- nbclust_args:

  Named list of arguments forwarded to
  [`NbClust::NbClust`](https://rdrr.io/pkg/NbClust/man/NbClust.html).
  Sensible defaults are set internally (`min.nc = 2`, `max.nc = 10`,
  `method = "ward.D2"`, `index = "all"`).

- thresh_metric:

  Character. Default `"sensitivity"`. The
  [`dismo::threshold`](https://rdrr.io/pkg/dismo/man/threshold.html)
  metric used to cut the future SDM into a binary surface.

## Value

A named list:

- `clusters_raster`:

  Smoothed `SpatRaster` of consensus cluster IDs.

- `clusters_sf`:

  `sf` polygons with column `ID`.

- `suitable_habitat`:

  `SpatRaster` - binary suitable habitat under future conditions.

- `novel_mask`:

  `SpatRaster` - `1` where MESS \< threshold.

- `mess`:

  `SpatRaster` - raw MESS scores.

- `stability`:

  `SpatRaster` - draw-based cluster-assignment stability under future
  climate. Fraction of posterior draws agreeing on each cell's modal
  cluster. Values in (0, 1); low values flag cells where beta
  uncertainty produces ambiguous cluster assignment in future space.
  Unlike the current-era stability surface, this is computed directly in
  future climate space and is not a transfer from
  [`PosteriorCluster()`](https://sagesteppe.github.io/safeHavens/reference/PosteriorCluster.md).

- `changes`:

  `data.frame` - per-cluster area and centroid-shift metrics.

- `novel_similarity`:

  `data.frame` - nearest existing cluster and mean silhouette width for
  each novel cluster. Zero-row if none.

## Details

1.  Rescale future predictors via `RescaleRasters_bayes` to match the
    feature space that
    [`PosteriorCluster()`](https://sagesteppe.github.io/safeHavens/reference/PosteriorCluster.md)
    used when training the consensus KNN.

2.  Run [`dismo::mess`](https://rdrr.io/pkg/dismo/man/mess.html) to
    identify novel climate cells.

3.  Predict known-climate cells with the consensus KNN classifier from
    [`PosteriorCluster()`](https://sagesteppe.github.io/safeHavens/reference/PosteriorCluster.md).

4.  If novel cells exist and `cluster_novel = TRUE`, cluster them
    independently with `NbClust` / KNN (`cluster_novel_areas`).

5.  Sample points from every cluster (existing + novel), build a single
    tree, and extract nearest-existing-cluster relationships via
    silhouette (`analyze_cluster_relationships`).

6.  Project future stability by replaying posterior beta draws on the
    future climate surface (`project_future_draws`). Stability here is a
    genuinely forward-looking quantity: the fraction of draws in which
    each cell's modal cluster assignment was agreed upon. This sidesteps
    the concern of transferring a current-era stability surface trained
    in current climate space onto a future surface that may differ
    substantially.

7.  Polygonise and calculate area / centroid changes.

## See also

[`projectClusters()`](https://sagesteppe.github.io/safeHavens/reference/projectClusters.md)
for the elastic-net equivalent,
[`PosteriorCluster()`](https://sagesteppe.github.io/safeHavens/reference/PosteriorCluster.md),
[`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md)
