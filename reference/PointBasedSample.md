# Generate a sampling grid based off of regularly sampled points across the species range.

This function utilizes a regular, or nearly so in the case of existing
collections, grid of points to develop a sampling scheme or n polygons.

## Usage

``` r
PointBasedSample(polygon, n = 20, collections, reps = 100, BS.reps = 9999)
```

## Arguments

- polygon:

  the input sf polygon, i.e. species range or administrative unit, where
  sampling is desired.

- n:

  Numeric. The total number of desired collections. Defaults to 20.

- collections:

  an sf point geometry data set of where existing collections have been
  made.

- reps:

  further arguments passed to np.boot

- BS.reps:

  number of bootstrap replicates for evaluating results.

## Value

A list containing two objects, the first the results of bootstrap
simulations. The second an sf dataframe containing the polygons with the
smallest amount of variance in size.

## Examples

``` r
#' Utilize a grid based stratified sample for drawing up polygons
ri <- spData::us_states |>
  dplyr::select(NAME) |>
  dplyr::filter(NAME == 'Rhode Island') |>
  sf::st_transform(32617)
  
 system.time(
  out <- PointBasedSample(polygon = ri, reps = 10, BS.reps = 10) # set very low for example
 )
#>    user  system elapsed 
#>   0.596   0.003   0.599 
# the function is actually very fast; 150 voronoi reps, with 9999 BS should only take about
# 2 seconds per species so not much concern on the speed end of things!
head(out$SummaryData)
#>                  Metric    Value
#> 1     variance.observed 16538241
#> 2        quantile.0.001 16541620
#> 3             lwr.95.CI 16538241
#> 4             upr.95.CI 17219398
#> 5    Voronoi.reps.asked       10
#> 6 Voronoi.reps.received        6
plot(out$Geometry)

```
