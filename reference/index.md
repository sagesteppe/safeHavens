# Package index

## Sampling

User facing functions for implementing the sampling schemes

- [`GridBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/GridBasedSample.md)
  : Create hexagonal grid based polygons over a geographic range
- [`PointBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/PointBasedSample.md)
  : Generate a sampling grid based off of regularly sampled points
  across the species range.
- [`EqualAreaSample()`](https://sagesteppe.github.io/safeHavens/reference/EqualAreaSample.md)
  : Create equal area polygons over a geographic range
- [`OpportunisticSample()`](https://sagesteppe.github.io/safeHavens/reference/OpportunisticSample.md)
  : Design additional collections around already existing collections
- [`IBDBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/IBDBasedSample.md)
  : Sample a species based on Isolation by Geographic Distance (IBD)
- [`EcoregionBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EcoregionBasedSample.md)
  : Determine the sample size for n ecoregions in an area
- [`EnvironmentalBasedSample()`](https://sagesteppe.github.io/safeHavens/reference/EnvironmentalBasedSample.md)
  : Create environmental and spatial clusters for targeting collection
  areas
- [`PrioritizeSample()`](https://sagesteppe.github.io/safeHavens/reference/PrioritizeSample.md)
  : Determine which areas of a sample unit should be prioritized

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

## Rare Plant Sampling

Optimization based site selection for known populations

- [`greatCircleDistance()`](https://sagesteppe.github.io/safeHavens/reference/greatCircleDistance.md)
  : Haversine Distance Calculation
- [`maximizeDispersion()`](https://sagesteppe.github.io/safeHavens/reference/maximizeDispersion.md)
  : Maximize Dispersion Site Selection

## Miscellaneous

Data prep functions

- [`TestGridSizes()`](https://sagesteppe.github.io/safeHavens/reference/TestGridSizes.md)
  : Get an estimate for how many grids to draw over a species range
