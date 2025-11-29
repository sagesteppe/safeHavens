# Rescale a raster stack to reflect the beta coefficients from a glmnet model

These rescaled rasters can then be used for clustering, and predicting
the results of cluster analysis back into space for a final product.

## Usage

``` r
RescaleRasters(model, predictors, training_data, pred_mat)
```

## Arguments

- model:

  the final output model from glmnet from `elasticSDM`

- predictors:

  the raster stack to use for the process from `elasticSDM`

- training_data:

  the same data that went into the glmnet model, this is used for
  calculating variance which is required for the scaling process. From
  `elasticSDM`

- pred_mat:

  the Prediction matrix from `elasticSDM`

## Value

A list with two objects. 1) The rescaled raster stack. 2) A table of
both standardized and unstandardized coefficients from the glmnet model.
