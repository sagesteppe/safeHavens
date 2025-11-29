# Get an estimate for how many grids to draw over a species range

This function uses the dimensions of a species grid to estimate how many
grids would need to be added in the x and y directions to cover it with
20 grid cells of roughly equal areas.

## Usage

``` r
TestGridSizes(target)
```

## Arguments

- target:

  a species range as a simple feature (`sf`) object.

## Value

A dataframe with testing results for each grid combination. A user needs
to select the optimal grid size based on a tradeoff with minimizing
variance, without creating too many grids which will need to be erased.
In the Rhode Island example I would use the 'Original' option which asks
for 4 x grids and 7 y grids.

## Examples

``` r
ri <- spData::us_states |>
dplyr::select(NAME) |>
   dplyr::filter(NAME == 'Rhode Island') |>
   sf::st_transform(32617)

sizeOptions <- TestGridSizes(ri)
head(sizeOptions)
#>       Name Grids    Variance GridNOx GridNOy
#> 1 Smallest    49    9.242133       6       9
#> 2  Smaller    37  355.799568       5       8
#> 3 Original    25 1038.567100       4       7
#> 4   Larger    16 1245.553777       3       6
#> 5  Largest    11 1322.092849       2       5
```
