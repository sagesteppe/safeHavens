# Fit a Bayesian spatial GLMM as a species distribution model

A Bayesian alternative to
[`elasticSDM()`](https://sagesteppe.github.io/safeHavens/reference/elasticSDM.md)
using `brms`. Fits a binomial GLMM with a Gaussian process spatial
random effect over the full pairwise distance matrix (via a kernel
covariance, no mesh approximation). All environmental predictors are
retained; shrinkage toward zero is handled by priors rather than
variable selection.

Users should either (a) remove highly correlated predictors beforehand,
or (b) pass PCA axes as predictors (see `pca_predictors` argument),
which is the recommended workflow.

The function returns an object structurally compatible with the
`safeHavens` downstream workflow: the same named list slots are
populated where meaningful, with Bayesian-specific additions (posterior
draws, LOO-CV).

## Usage

``` r
bayesianSDM(
  x,
  predictors,
  planar_projection,
  quantile_v = 0.025,
  prior_type = c("horseshoe", "normal", "student"),
  prior_scale = 1,
  pca_predictors = TRUE,
  pca_axes = 5,
  gp_scale_prior = NULL,
  resample = FALSE,
  feature_selection = c("ffs", "none"),
  vif = TRUE,
  min_ffs_var = 5,
  chains = 4,
  iter = 5000,
  warmup = 2000,
  cores = 4,
  k = 5,
  seed = 42,
  backend = "cmdstanr",
  fact = 3,
  ...
)
```

## Arguments

- x:

  An `sf` object of occurrence points. Must have an `occurrence` column
  coded 0/1, OR the function will assume all rows are presences and
  generate pseudo-absences (matching
  [`elasticSDM()`](https://sagesteppe.github.io/safeHavens/reference/elasticSDM.md)
  behaviour).

- predictors:

  A `terra` SpatRaster stack of environmental predictors.

- planar_projection:

  Numeric EPSG code or proj4 string for a planar (metre-unit) CRS. Used
  for spatial thinning and to derive GP coordinates. See
  [`elasticSDM()`](https://sagesteppe.github.io/safeHavens/reference/elasticSDM.md)
  for guidance.

- quantile_v:

  Numeric thinning quantile passed to
  [`spThin::thin()`](https://rdrr.io/pkg/spThin/man/thin.html). Defaults
  to `0.025`. Set to `0.001` for essentially no thinning.

- prior_type:

  Character. Prior family placed on the fixed-effect environmental
  coefficients. One of:

  `"horseshoe"`

  :   Regularised horseshoe (Piironen & Vehtari 2017). Best when many
      predictors may be irrelevant. Requires the expected number of
      non-zero coefficients via `p0` (see `...`).

  `"normal"`

  :   `Normal(0, prior_scale)`. Simple and fast; reasonable when
      predictors are pre-selected or are PCA axes.

  `"student"`

  :   Student-t(3, 0, prior_scale). Heavier tails than Normal; a middle
      ground.

  Defaults to `"horseshoe"`.

- prior_scale:

  Numeric. Scale parameter for `"normal"` and `"student"` priors.
  Ignored for `"horseshoe"`. Defaults to `1` (weakly informative on the
  log-odds scale when predictors are standardised).

- pca_predictors:

  Logical. If `TRUE`, replace the raw predictor stack with its principal
  components before fitting. The number of axes retained is controlled
  by `pca_axes`. Defaults to `FALSE`.

- pca_axes:

  Integer. Number of PCA axes to retain when `pca_predictors = TRUE`.
  Defaults to `5`.

- gp_scale_prior:

  A `brms` prior object for the GP length-scale parameter, or `NULL` to
  use brms default. The default `NULL` uses `inv_gamma(3, 1)`, which is
  weakly informative and appropriate for most ecological datasets.
  Override only if you have strong prior knowledge about the spatial
  range of your species (e.g.,
  `prior(inv_gamma(5, 2), class = lscale, coef = "")` for tighter
  spatial autocorrelation).

- resample:

  Boolean, Defaults to FALSE. Used to place 15% of the requested points
  in areas undersampled by sdm::background functions.

- feature_selection:

  Character. Variable selection method to apply before Bayesian model
  fitting. One of:

  `"ffs"`

  :   Forward feature selection via CAST::ffs(). Uses spatial CV folds
      to select variables. Fast and spatially-aware. Default.

  `"none"`

  :   No feature selection; use all predictors (or all PCA axes). Relies
      on horseshoe prior for shrinkage.

- vif:

  Boolean, default TRUE, Whether to run usdm::vifcor to remove highly
  collinear features, default theta used.

- min_ffs_var:

  Integer. Minium number of ffs vars to start with.

- chains:

  Integer. Number of MCMC chains. Defaults to `4`.

- iter:

  Integer. Total iterations per chain (including warmup). Defaults to
  `5000`.

- warmup:

  Integer. Warmup iterations per chain. Defaults to `1000`.

- cores:

  Integer. Parallel cores. Defaults to
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html).

- k:

  Integer. Number of spatial CV folds (CAST `knndm`). Defaults to `5`.

- seed:

  Integer. Random seed for reproducibility. Defaults to `42`.

- backend:

  Character. Stan backend passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html). One
  of `"cmdstanr"` (recommended, faster) or `"rstan"`. Defaults to
  `"cmdstanr"`.

- fact:

  Numeric, default 2.0. Factor to multiple the number of occurrence
  records by to generate the number of background (absence) points. \#'

- ...:

  Additional arguments forwarded to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
  Also accepts `p0` (integer, expected non-zero predictors) for the
  horseshoe prior; defaults to `floor(ncol(pred_matrix) / 2)`.

## Value

A named list with the following elements, mirroring
[`elasticSDM()`](https://sagesteppe.github.io/safeHavens/reference/elasticSDM.md):

- `RasterPredictions`:

  `SpatRaster` of posterior mean predicted occurrence probability.

- `RasterPredictions_sd`:

  `SpatRaster` of posterior SD (uncertainty).

- `Predictors`:

  The predictor raster stack actually used (may be PCA axes if
  `pca_predictors = TRUE`).

- `PCNM`:

  Always `NULL`; slot retained for workflow compatibility.

- `Model`:

  The fitted `brmsfit` object.

- `CVStructure`:

  CAST `knndm` object with fold indices.

- `LOO`:

  `loo` object from approximate leave-one-out cross-validation via
  PSIS-LOO. Replaces `ConfusionMatrix`.

- `ConfusionMatrix`:

  A simple threshold-based confusion matrix (threshold = 0.5 on
  posterior mean) on the held-out test fold, for comparability with
  [`elasticSDM()`](https://sagesteppe.github.io/safeHavens/reference/elasticSDM.md)
  output.

- `TrainData`:

  `sf` object used for model fitting.

- `TestData`:

  `sf` object held out for evaluation.

- `PredictMatrix`:

  The numeric matrix of predictor values used for the spatial raster
  prediction.

- `PCAModel`:

  The `prcomp` object if `pca_predictors = TRUE`, otherwise `NULL`.
  Retained so users can project new data onto the same axes.

- `Diagnostics`:

  Named list: `Rhat` (max R-hat across all parameters), `BulkESS` (min
  bulk ESS), `TailESS` (min tail ESS). Flags convergence problems
  automatically.

## See also

[`elasticSDM()`](https://sagesteppe.github.io/safeHavens/reference/elasticSDM.md)
for the elastic-net alternative,
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html),
[`CAST::knndm()`](https://hannameyer.github.io/CAST/reference/knndm.html)

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
library(terra)

x <- read.csv(file.path(system.file(package = "dismo"), "ex", "bradypus.csv"))
x <- sf::st_as_sf(x[, c("lon", "lat")], coords = c("lon", "lat"), crs = 4326)

files <- list.files(
  path    = file.path(system.file(package = "dismo"), "ex"),
  pattern = "grd",
  full.names = TRUE
)
predictors <- terra::rast(files)

# Recommended workflow: PCA axes as predictors, horseshoe prior
result <- bayesianSDM(
  x                 = x,
  predictors        = predictors,
  planar_projection = 5070,
  pca_predictors    = TRUE,
  pca_axes          = 5,
  prior_type        = "horseshoe",
  chains            = 4,
  iter              = 2000
)

terra::plot(result$RasterPredictions)
terra::plot(result$RasterPredictions_sd)  # uncertainty surface
print(result$Diagnostics)
} # }
```
