# Modify the output rasters from the SDM process to better match target goals

This is the last 'analytical' portion of the SD modelling process. It
will produce binary (Yes/ NA) rasters of species suitable habitat based
on a three step process. The first step uses dismo::thresholds to
determine what feature of the raster we want to maximuize, in this case
I want a raster to be *more* than less likely to capture a presence than
omit it. Then we can subset the predicted habitat, which has not been
dispersed too by comparing nearest neighbor distances of the observed
points. Finally, we can add back in areas to the raster where we know
that the species has been observed, hopefully these are not missed in
the original SDM, but there are always some suspicious points which are
difficult to fit in light of the inertia of the rest of the species.

## Usage

``` r
PostProcessSDM(
  rast_cont,
  test,
  train,
  thresh_metric = "sensitivity",
  quant_amt = 0.25,
  planar_projection
)
```

## Arguments

- rast_cont:

  the raw unaltered (except masked) raster predictions
  (x\$RasterPredictions).

- test:

  the test data partition from the `elasticSDM` function (x\$TestData).

- train:

  the train data partition from the `elasticSDM` function
  (x\$TrainData).

- thresh_metric:

  ?dismo::threshold for all options, defaults to 'sensitivity'

- quant_amt:

  the quantile of nearest neighbors distance to use for steps 2 and 3.
  defaults to 0.25, using the median nearest neighbor distance of 10
  bootstrapping replicates for estimating a buffer to restrict the SDM
  surface too, and the minimum of the 10 bootstrap reps for adding
  surface to presence points which were not placed in binary suitable
  habitat.

- planar_projection:

  Numeric, or character vector. An EPSG code, or a proj4 string, for a
  planar coordinate projection, in meters, for use with the function.
  For species with very narrow ranges a UTM zone may be best (e.g. 32611
  for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). Otherwise
  a continental scale projection like 5070 See
  https://projectionwizard.org/ for more information on CRS. The value
  is simply passed to sf::st_transform if you need to experiment.

## Value

A list containing two options. 1) A spatraster with 4 layers, A) the
continuous probabilities of suitable habitat feed in from `elasticSDM`,
B) this raster in binary format based on the specified thresholding
statistic, C) the binary raster from B + with habitat clipped to the
buffer distances determined by measuring nearest neighbor distances and
thresholding at a quantile D) the binary raster from C, adding the same
distance to points which were initially in cells classified by
thresholding as not having suitable habitat. D' is the general basis for
all future steps, but either B, or C serve as alternatives. 2) All
threshold statistics calculated by dismo as a dataframe.
