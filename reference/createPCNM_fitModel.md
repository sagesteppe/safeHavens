# Create global/regional PCNM surfaces and fit elastic net regression to all covariates

This function does most of the lifting in the SDM workflow. It will
create PCNM/MEM surfaces and subset them to the local/global maps. It
can then use Thin plate regression to predict those onto an actual
raster surface which can be used for prediction downstream.

## Usage

``` r
createPCNM_fitModel(x, planar_proj, ctrl, indices_knndm, sub, test)
```

## Arguments

- x:

  should be the training data as an sf/tibble/dataframe

- planar_proj:

  Numeric, or character vector. An EPSG code, or a proj4 string, for a
  planar coordinate projection, in meters, for use with the function.
  For species with very narrow ranges a UTM zone may be best (e.g. 32611
  for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). Otherwise
  a continental scale projection like 5070 See
  https://projectionwizard.org/ for more information on CRS. The value
  is simply passed to sf::st_transform if you need to experiment.

- ctrl:

  the control object created by character in the SDM function.

- indices_knndm:

  from sdm function

- sub:

  the subset predictors from elasticSDM

- test:

  the test data partition from elasticSDM

## Details

It will then use cross validation to determine a suitable glmnet model
alpha and lambda, and fit them using glmnet. It returns three objects
which are spit out into the environment, 1) pcnm, surfaces for only
those eigenvectors used in the glmnet model (including if shrunk out),
2) the glmnet model 3) all fitting information from carets process.
