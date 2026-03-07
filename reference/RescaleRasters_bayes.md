# Rescale a raster stack to reflect posterior beta coefficients from a brms SDM

Bayesian analogue of
[`RescaleRasters()`](https://sagesteppe.github.io/safeHavens/reference/RescaleRasters.md)
for use with
[`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md)
output. Each predictor raster is centred and scaled using training-data
moments, then multiplied by the absolute value of its posterior beta
coefficient, so that the resulting layers are on a common scale weighted
by ecological importance.

Because `brms` uses `autoscale = TRUE` by default for regularising
priors (horseshoe, normal, student), the posterior betas are already on
a standardised scale and are directly comparable across predictors - no
additional correction is needed beyond what
[`RescaleRasters()`](https://sagesteppe.github.io/safeHavens/reference/RescaleRasters.md)
applies to glmnet output.

## Usage

``` r
RescaleRasters_bayes(
  model,
  predictors,
  training_data,
  pred_mat,
  beta_summary = c("mean", "median", "Q2.5", "Q97.5"),
  include_uncertainty = FALSE,
  uncertainty_wt = 1
)
```

## Arguments

- model:

  A `brmsfit` object from
  [`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md).

- predictors:

  A `SpatRaster` stack. Should be the `$Predictors` element from
  [`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md)
  output (i.e. already PCA-rotated if `pca_predictors` was `TRUE`).

- training_data:

  An `sf` object. The `$TrainData` element from
  [`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md)
  output. Used to compute per-variable mean and SD for centring/scaling
  the raster.

- pred_mat:

  A data frame or matrix. The `$PredictMatrix` element from
  [`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md)
  output. Column names must match environmental predictor names (GP
  coordinate columns `gp_x`/`gp_y` are ignored automatically).

- beta_summary:

  Character. Which posterior summary to use as the beta weight. One of
  `"mean"`, `"median"`, `"Q2.5"`, `"Q97.5"`. Defaults to `"mean"`.

- include_uncertainty:

  Logical. Whether to append a layer encoding posterior uncertainty of
  the linear predictor. Defaults to `FALSE`.

- uncertainty_wt:

  Numeric. Weight applied to the uncertainty layer relative to the
  maximum environmental layer range, analogous to `coord_wt` in
  [`EnvironmentalBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EnvironmentalBasedSample.md).
  Only used when `include_uncertainty = TRUE`. Defaults to `1`.

## Value

A list with three elements:

- `RescaledPredictors`:

  `SpatRaster` of rescaled, beta-weighted predictor layers. Pass this as
  `pred_rescale` to
  [`EnvironmentalBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EnvironmentalBasedSample.md).

- `BetaCoefficients`:

  Data frame of posterior summaries for all fixed effects (columns:
  `Variable`, `Estimate` (mean), `Est.Error` (SD), `Q2.5`, `Q97.5`,
  `BetaWeight` (the value actually used for scaling)). Intercept and GP
  parameters are excluded.

- `UncertaintyLayer`:

  `SpatRaster` of propagated posterior SD, or `NULL` if
  `include_uncertainty = FALSE`.

## Posterior summary choices

The `beta_summary` argument controls which posterior summary statistic
is used as the point-estimate beta weight:

- `"mean"`:

  Posterior mean. Minimum MSE estimator; default and recommended for
  most use cases.

- `"median"`:

  Posterior median. More robust when horseshoe priors produce
  heavy-tailed marginals for near-zero coefficients.

- `"Q2.5"` / `"Q97.5"`:

  Credible interval bounds. Useful for sensitivity analysis (e.g.
  conservative lower-bound weighting).

## Uncertainty layer

When `include_uncertainty = TRUE` an additional raster layer is appended
whose values are the *posterior SD of the linear predictor* at each
cell, propagated through all environmental betas. This encodes where the
model is most uncertain about environmental conditions,

## See also

[`RescaleRasters()`](https://sagesteppe.github.io/safeHavens/reference/RescaleRasters.md)
for the glmnet equivalent,
[`bayesianSDM()`](https://sagesteppe.github.io/safeHavens/reference/bayesianSDM.md),
[`EnvironmentalBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EnvironmentalBasedSample.md)
