# Sample spatial zones within a species range

Intersect a vector data file of spatial zones (ecoregions, provisional
seed transfer zones, or other spatial partitions) with the range of a
focal taxon and select `n` zones to sample. If fewer than n zones exist,
extra samples are allocated using the specified method. If more than n
zones exist, zones are selected using the specified method.

## Usage

``` r
PolygonBasedSample(
  x,
  zones,
  zone_key,
  n = 20,
  decrease_method = c("Largest", "Smallest", "Most", "Assist-warm", "Assist-drier"),
  increase_method = c("Largest", "Smallest", "Most", "Assist-warm", "Assist-drier"),
  warmest_col = NULL,
  precip_col = NULL
)
```

## Arguments

- x:

  sf object. Species range as a simple feature.

- zones:

  sf object. Spatial zones vector data (ecoregions, PSTZs, etc.).

- zone_key:

  Character. Column name identifying unique zones (required).

- n:

  Numeric. Desired total number of samples. Default = 20.

- decrease_method:

  Character. Method when n \< number of zones. One of: "Largest",
  "Smallest", "Most", "Assist-warm", "Assist-drier". Default =
  "Largest".

- increase_method:

  Character. Method when n \> number of zones. One of: "Largest",
  "Smallest", "Most", "Assist-warm", "Assist-drier". Default =
  "Largest".

- warmest_col:

  Character. Column name for warmest temperature metric (required for
  "Assist-warm" method).

- precip_col:

  Character. Column name for precipitation metric (required for
  "Assist-drier" method).

## Value

sf object with selected zones/polygons and an `allocation` column
indicating number of samples per polygon.

## Details

Simple features can store polygon data as 'MULTIPOLYGON' (all polygons
of a class stored collectively) or 'POLYGON' (each individual polygon is
a unique entry). This function will cast MULTIPOLYGONS to POLYGONS as
needed, but pre-casting will improve performance.

Available methods for zone selection:

- **Largest**: Select zones by total area (descending)

- **Smallest**: Select zones by total area (ascending)

- **Most**: Select zones with most polygons (highest fragmentation)

- **Assist-warm**: Select warmest zones (requires `warmest_col`)

- **Assist-drier**: Select driest zones (requires `precip_col`)

## Examples

``` r
if (FALSE) { # \dontrun{
library(tidyverse)

sr_mat <- rbind(
  c(0,0), c(10,0), c(10,10), c(0,10), c(0,0)
)
sr_poly <- sf::st_polygon(list(sr_mat))
x <- sf::st_sf(id = 1, geometry = sf::st_sfc(sr_poly))

rm(sr_mat, sr_poly)

zone_polys = data.frame(
  ## randomly generate some points in XY space. 
  x = runif(10, min = -2, max = 12),
  y = runif(10, min = -2, max = 12)
) |>
  # conver to spatial points
  sf::st_as_sf(coords = c('x', 'y')) |>
  ## allocate XY space to it's nearest point
  sf::st_union() |>
  sf::st_voronoi() |>
  ## extract the contiguous pieces of XY space around points
  sf::st_collection_extract('POLYGON') |>
  sf::st_as_sf() |>
  ## make up seed zones on the fly, assign multiple polygons to some zones. 
  dplyr::mutate(pstz_key = sample(LETTERS[1:7], size = 10, replace = T)) |>
  dplyr::rename('geometry' = x) |>
  sf::st_crop(x)

bp <- ggplot2::ggplot(x) + 
  ggplot2::geom_sf(fill = NA, lwd = 2) + 
  ggplot2::geom_sf(data = zone_polys, ggplot2::aes(fill = pstz_key)) 

bp + 
  ggplot2::geom_sf_label(data = zone_polys, ggplot2::aes(label = pstz_key))

###################################################################### 
# example #1: request same numer of samples as zones - all zones returned. 
 res1 <- PolygonBasedSample(
   x = x, 
   n = length(unique(zone_polys[['pstz_key']])), 
   zones = zone_polys, 
   zone_key  = "pstz_key",
   increase_method = "Most"
 )

bp +
  geom_sf(data = res1, alpha = 0.9) + 
  geom_sf_label(data = pstz,aes(label = pstz_key)) 

## note that we get the largest polygon from EACH group to sample from. 

##################################################################### 
# Example #2: request fewer samples than zones -> subset by method - choosing largest by area 

res2 <- PolygonBasedSample(
   x = x, n = 3, zones = zone_polys, zone_key  = "pstz_key", increase_method = "Largest"
)

res2 |>
  group_by(pstz_key) |>
  mutate(total_area = sum(poly_area)) |>
  sf::st_drop_geometry() |>
  arrange(-total_area)|>
  knitr::kable()

bp + # picks, the three largest 
  geom_sf(data = res2, alpha = 0.9) + 
  geom_sf_label(data = zone_polys, aes(label = pstz_key))

####################################################################### 
# Example #3: request fewer samples than zones -> subset by method - choosing smallest by area 
res3 <- PolygonBasedSample(
  x = x, n = 3, zones = zone_polys, zone_key = "pstz_key", increase_method = "Smallest")

res3 |>
  group_by(pstz_key) |>
  mutate(total_area = sum(poly_area)) |>
  sf::st_drop_geometry() |>
  arrange(total_area) |>
  knitr::kable()

## returns the largest polygon (poly_area) within the `pstz_key` group, ranked by (total_area)

bp + # picks, the n smallest - too small to see sometimes
  geom_sf(data = filter(res3, allocation == 0), alpha = 0.9) + 
  geom_sf_label(data = zone_polys, aes(label = pstz_key)) 

####################################################################
# Example #4: request more samples than zones -> allocate extras to Largest pSTZs

## note that is really a rounding rule - 'Largest' favors giving extra collections to the largest 
## polygons while 'smallest' favors giving them smaller polygons. It is really mostly for edge cases
## and the two will generally behave similarly on contrived examples. 
res4 <- PolygonBasedSample(
   x = x, n = 12, zones = zone_polys, zone_key = "pstz_key", increase_method = "Largest")
  
res4 |>
  group_by(pstz_key) |>
  summarize(total_area = sum(poly_area),  Total_Allocation = sum(allocation)) |>
  sf::st_drop_geometry() |>
  arrange(-total_area) |>
  knitr::kable()

bp + 
  theme(legend.position = 'none') + 
  geom_sf(data = res4, aes(fill = as.factor(allocation))) + 
  geom_sf_label(data = res4, aes(label = allocation)) 

####################################################################
# Example #5: request more samples than zones -> allocate extras to pSTZs with most polygons
res5 <- PolygonBasedSample(
   x = x, n = 14, zones = zone_polys, zone_key = "pstz_key", increase_method = "Most")
  
res5 |>
  group_by(pstz_key) |>
  summarize(Count = n(), Total_Allocation = sum(allocation)) |>
  sf::st_drop_geometry() |>
  arrange(-Count) |>
  knitr::kable()
} # }
```
