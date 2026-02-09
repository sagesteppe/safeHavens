# Rescale future climate predictors using current model coefficients

Standardises future predictors using the mean and SD from *current*
climate, then applies the glmnet beta weights. This is the same
weighting applied upstream in `RescaleRasters` — here we just do it
against a different raster stack.

## Usage

``` r
rescaleFuture(
  model,
  future_predictors,
  current_predictors,
  training_data,
  pred_mat
)
```

## Arguments

- model:

  glmnet model object from `elasticSDM_noPCNM()$Model`.

- future_predictors:

  SpatRaster of future climate variables. Names must match those
  retained in `model`.

- current_predictors:

  SpatRaster of current climate (provides the standardisation
  parameters).

- training_data:

  the same data that went into the glmnet model, this is used for
  calculating variance which is required for the scaling process. From
  `elasticSDM`

- pred_mat:

  the Prediction matrix from `elasticSDM`

## Value

SpatRaster with rescaled future predictors (one layer per non-zero
coefficient).
