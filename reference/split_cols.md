# split and extract the temperature values from Tmin and AHM columns

Programmed for the Bower provisional seed zone products, a helper
function for separating and recovering the values from the columns in
the data.

## Usage

``` r
split_cols(dat, y, sep = "-")
```

## Arguments

- dat:

  data frame with the columns required to split.

- y:

  character. column name to split.

- sep:

  character. separator between the values to split on. Default is '-'.

## Examples

``` r
df = data.frame(
  'Tmin_class' = c('10 - 15 Deg. F.', '15 - 20 Deg. F.', '> 55 Deg. F.' ),
  'AHM_class' = c('2 - 3', '6 - 12', '3 - 6')
)
split_cols(df, 'Tmin_class')
#>   lower upper median range
#> 1    10    15   12.5     5
#> 2    15    20   17.5     5
#> 3    55    55   55.0     0
split_cols(df, 'AHM_class')
#>   lower upper median range
#> 1     2     3    2.5     1
#> 2     6    12    9.0     6
#> 3     3     6    4.5     3
```
