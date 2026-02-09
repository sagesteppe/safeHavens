# Package index

## Sampling

User facing functions for implementing the sampling schemes

- [`PointBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/PointBasedSample.md)
  : Generate a sampling grid based off of regularly sampled points
  across the species range.
- [`EqualAreaSample()`](https://sagesteppe.github.io/safeHavens/reference/EqualAreaSample.md)
  : Create equal area polygons over a geographic range
- [`OpportunisticSample()`](https://sagesteppe.github.io/safeHavens/reference/OpportunisticSample.md)
  : Design additional collections around already existing collections
- [`KMedoidsBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/KMedoidsBasedSample.md)
  : K Medoids Based Sample Site Selection Select a subset of sites that
  maximize spatial dispersion of sites using k-medioids clustering.
- [`IBDBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/IBDBasedSample.md)
  : Sample a species based on Isolation by Geographic Distance (IBD)
- [`PolygonBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/PolygonBasedSample.md)
  : Sample spatial zones within a species range
- [`EnvironmentalBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EnvironmentalBasedSample.md)
  : Create environmental and spatial clusters for targeting collection
  areas
- [`PrioritizeSample()`](https://sagesteppe.github.io/safeHavens/reference/PrioritizeSample.md)
  : Determine which areas in a sample unit should be prioritized

## Isolation by Resistance Sampling Prep

Preparation for Isolation by Resistance Sampling

- [`buildResistanceSurface()`](https://sagesteppe.github.io/safeHavens/reference/buildResistanceSurface.md)
  : Create a simple, theoretical, raster surface modelling Isolation by
  Distance.
- [`populationResistance()`](https://sagesteppe.github.io/safeHavens/reference/populationResistance.md)
  : Identify clusters of populations least separated by landscape
  resistance
- [`IBRSurface()`](https://sagesteppe.github.io/safeHavens/reference/IBRSurface.md)
  : Sample a species based on Isolation by Resistance Distance (IBR)

## Environmental Sampling Prep

Preparation for Environmental Based Sampling

- [`elasticSDM()`](https://sagesteppe.github.io/safeHavens/reference/elasticSDM.md)
  : Create a quick SDM using elastic net regression

- [`PostProcessSDM()`](https://sagesteppe.github.io/safeHavens/reference/PostProcessSDM.md)
  : Modify the output rasters from the SDM process to better match
  target goals

- [`RescaleRasters()`](https://sagesteppe.github.io/safeHavens/reference/RescaleRasters.md)
  : Rescale a raster stack to reflect the beta coefficients from a
  glmnet model

- [`writeSDMresults()`](https://sagesteppe.github.io/safeHavens/reference/writeSDMresults.md)
  :

  Save the results of SDMs from the `elasticSDM`, `RescaleRasters` and
  `PostProcessSDM` function in safeHavens

## Predictive Provenance Prep

Preparation for Predictive Provenance

- [`rescaleFuture()`](https://sagesteppe.github.io/safeHavens/reference/rescaleFuture.md)
  : Rescale future climate predictors using current model coefficients
- [`projectClusters()`](https://sagesteppe.github.io/safeHavens/reference/projectClusters.md)
  : Project current environmental clusters onto a future climate
  scenario

## Rare Plant Sampling

Optimization based site selection helpers for known populations

- [`greatCircleDistance()`](https://sagesteppe.github.io/safeHavens/reference/greatCircleDistance.md)
  : Haversine Distance Calculation
- [`split_cols()`](https://sagesteppe.github.io/safeHavens/reference/split_cols.md)
  : split and extract the temperature values from Tmin and AHM columns
