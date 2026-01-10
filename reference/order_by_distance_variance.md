# Order zones by minimizing distance variance

Order a set of spatial zones based on minimizing distance variance

## Usage

``` r
order_by_distance_variance(x, metric = c("var", "sd", "energy", "cv"))
```

## Arguments

- x:

  an sf object containing the zones to be ordered

- metric:

  character. The metric to minimize when ordering zones. Options are
  "var" (variance), "sd" (standard deviation), "energy" (sum of squared
  distances), and "cv" (coefficient of variation).

## Value

A numeric vector representing the order of zones based on the specified
distance variance metric
