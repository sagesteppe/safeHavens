# Cluster environmental space using posterior beta draws

Implements option D of the Bayesian clustering workflow: draws K sets of
beta coefficients from the posterior, produces a weighted environmental
raster for each draw, samples points from each, and records how often
pairs of sample points are assigned to the same cluster. The resulting
co-occurrence probability matrix is then used to derive a consensus
clustering.

Cells that co-cluster with probability near 1 across all draws are
robustly similar in environmental space. Cells near boundaries that
shift across draws — where the model is uncertain which variables matter
most — will have intermediate co-occurrence probability. This
boundary-instability surface is returned separately as `StabilityRaster`
and is directly meaningful for targeting collection effort.

## Usage

``` r
PosteriorCluster(
  model,
  predictors,
  f_rasts,
  pred_mat,
  training_data,
  n_draws = 100,
  n = 10,
  n_pts = 500,
  lyr = "occurrence_prob_mean",
  planar_proj,
  coord_wt = 2.5,
  consensus_method = c("hierarchical", "pam"),
  beta_draws = NULL,
  seed = 42
)
```

## Arguments

- model:

  A `brmsfit` object from
  [`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md).

- predictors:

  A `SpatRaster` stack (`$Predictors` from
  [`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md)).

- f_rasts:

  The rasters output from
  [`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md),
  used to define the spatial mask and domain. Specifically uses
  `$RasterPredictions`.

- pred_mat:

  Data frame or matrix (`$PredictMatrix` from
  [`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md)).

- training_data:

  `sf` object (`$TrainData` from
  [`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md)).

- n_draws:

  Integer. Number of posterior beta draws to cluster over. Defaults to
  `100`. Values between 50–200 give stable co-occurrence matrices for
  typical SDM datasets (\< 2000 obs, \< 20 vars).

- n:

  Integer. Number of clusters per draw. Passed to
  [`EnvironmentalBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EnvironmentalBasedSample.md)
  internals.

- n_pts:

  Integer. Points to sample per draw for clustering. Defaults to `500`.
  The same points are used across all draws (fixed spatial frame) so
  co-occurrence is directly comparable.

- lyr:

  Character. Layer name in `f_rasts` to use as the spatial domain mask.
  Defaults to `"Threshold"` (the
  [`PostProcessSDM()`](https://sagesteppe.github.io/safeHavens/reference/PostProcessSDM.md)
  output name).

- planar_proj:

  Numeric EPSG or proj4 string for planar projection.

- coord_wt:

  Numeric. Coordinate weight, as in
  [`EnvironmentalBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EnvironmentalBasedSample.md).
  Defaults to `2.5`.

- consensus_method:

  Character. How to derive the final clustering from the co-occurrence
  matrix. One of:

  `"hierarchical"`

  :   Average-linkage hierarchical clustering on `(1 - co_occurrence)`
      as a dissimilarity. Stable and interpretable. Default.

  `"pam"`

  :   Partitioning Around Medoids on the dissimilarity matrix. More
      robust to outliers than hierarchical.

- beta_draws:

  Optional matrix of pre-drawn posterior samples (rows = draws, cols =
  parameters, as from
  [`brms::as_draws_matrix()`](https://paulbuerkner.com/brms/reference/draws-brms.html)).
  If `NULL`, draws are taken internally. Supply this to reuse draws
  across multiple calls.

- seed:

  Integer. Random seed. Defaults to `42`.

## Value

A list:

- `Geometry`:

  `sf` polygon layer of consensus clusters (compatible with
  [`EnvironmentalBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EnvironmentalBasedSample.md)
  `$Geometry` output).

- `CoOccurrenceMatrix`:

  Numeric matrix (n_pts × n_pts). Each tupple is the proportion of draws
  in which points i and j were assigned to the same cluster.

- `StabilityRaster`:

  `SpatRaster` of per-cell cluster-assignment stability, predicted via
  KNN regression from the sample point stability scores. Values in
  (0,1); low values flag environmentally ambiguous boundary regions.

- `ConsensusRaster`:

  `SpatRaster` of final consensus (rank1) cluster IDs, geographically
  reordered.

- `Rank2Raster`:

  `SpatRaster` of second-most-frequent cluster assignments across draws.
  Where a point was always assigned to one cluster, rank2 = rank1.

- `Rank3Raster`:

  `SpatRaster` of third-most-frequent cluster assignments. Where fewer
  than 3 unique clusters were assigned, rank3 = rank1.

- `Top3Lookup`:

  Data frame with columns: `point_id`, `x`, `y`, `rank1_cluster`,
  `rank1_pct`, `rank2_cluster`, `rank2_pct`, `rank3_cluster`,
  `rank3_pct`, `uncertainty` (= 100 - rank1_pct). Tabular summary of
  cluster assignment frequencies at sample points.

- `DrawClusterings`:

  Integer matrix (n_pts × n_draws) of per-draw cluster assignments.
  Retained for diagnostics.

- `SamplePoints`:

  `sf` object of the fixed sample points used.

- `KNN_Cluster`:

  Trained KNN model for rank1 cluster assignment. Can be reused to
  classify new points.

- `KNN_Rank2`:

  Trained KNN model for rank2 clusters.

- `KNN_Rank3`:

  Trained KNN model for rank3 clusters.

- `KNN_Stability`:

  Trained KNN regression model for stability scores. Can be reused to
  predict stability at new locations.

## See also

[`RescaleRasters_bayes()`](https://sagesteppe.github.io/safeHavens/reference/RescaleRasters_bayes.md),
[`EnvironmentalBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EnvironmentalBasedSample.md),
[`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md)
