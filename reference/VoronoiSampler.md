# Make a voronoi sample of an area n times

Split an area up into n polygons of roughly equal area, optionally
removing some of the default points and replacing them with existing
collections to build the future collections around.

## Usage

``` r
VoronoiSampler(polygon, n, collections, reps)
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

  Numeric. The number of times to rerun the voronoi algorithm, the set
  of polygons with the most similar sizes, as measured using their
  variance of areas will be selected. Defaults to 150, which may
  accomplish around 100 succesful iterations.
