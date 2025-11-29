# Save the results of SDMs from the `elasticSDM`, `RescaleRasters` and `PostProcessSDM` function in safeHavens

This function is used to write out a wide range of values from the
`fitPredictOperationalize` process. It will create multiple
subdirectories within a user specified path. These include: 'Rasters'
where a raster stack of the four final rasters will go, 'Fitting' where
the details of model fitting from caret will be placed, 'Models' where
the final fit model will go, 'Evaluation' where all evaluation
statistics will be placed, 'Threshold' where results form
dismo::threshold will be placed.

## Usage

``` r
writeSDMresults(
  path,
  taxon,
  cv_model,
  pcnm,
  model,
  cm,
  coef_tab,
  f_rasts,
  thresh
)
```

## Arguments

- path:

  a root path where each of 5 folders will be created, if they do not
  exist.

- taxon:

  the name of the taxonomic entity for which the models were created.

- cv_model:

  the cross validation data from `elasticSDM`

- pcnm:

  the pcnm/mem rasters from `elasticSDM`

- model:

  the final glmnet model from `elasticSDM`

- cm:

  the confusion matrix from `elasticSDM`

- coef_tab:

  the coefficient table from `RescaleRasters`

- f_rasts:

  the final rasters from `RescaleRasters`

- thresh:

  threshold statistics from `PostProcessSDM`

## Value

all the above objects, or all objects above specified, are written to
disk.
