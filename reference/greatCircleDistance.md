# Haversine Distance Calculation

Calculate geographic distances on the geoid. This results in more
accurate distance calculations than a planar system. Function is mostly
used internally by `maximize_dispersion`

## Usage

``` r
greatCircleDistance(lat1, lon1, lat2, lon2)
```

## Arguments

- lat1:

  Double. column holding coords of 'focal' population

- lon1:

  Double. column holding coords of 'focal' population

- lat2:

  Double. column holding coords of 'non-focal' population

- lon2:

  Double. column holding coords of 'non-focal' population

## Details

calculate distances between sites (Haversine formula)

## Examples

``` r
n_sites <- 5 # number of known populations
 df <- data.frame(
   site_id = seq_len(n_sites),
   lat = runif(n_sites, 25, 30),
   lon = runif(n_sites, -125, -120)
 )

dist_mat <- sapply(1:nrow(df), function(i) {
   greatCircleDistance(
     df$lat[i], df$lon[i],
     df$lat, df$lon
   )
 })
#head(dist_mat)
```
