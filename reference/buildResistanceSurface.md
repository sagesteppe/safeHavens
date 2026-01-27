# Create a simple, theoretical, raster surface modelling Isolation by Distance.

Used in conjunction, preferably before `populationResistance` in the
`IBRBasedSample` workflow.

## Usage

``` r
buildResistanceSurface(
  base_raster,
  resistance_surface = NULL,
  oceans = NULL,
  lakes = NULL,
  rivers = NULL,
  tri = NULL,
  habitat = NULL,
  w_ocean = 1000,
  w_lakes = 200,
  w_rivers = 20,
  w_tri = 1,
  w_habitat = 1
)
```

## Arguments

- base_raster:

  SpatRaster. Base raster for the study area. Provides template geometry
  and resolution.

- resistance_surface:

  SpatRaster. Optional pre-computed resistance raster. If provided, the
  raster-building arguments are ignored.

- oceans:

  SpatRaster. Binary (0/1) raster for ocean cells. Used to increase
  movement cost.

- lakes:

  SpatRaster. Binary (0/1) raster for lakes.

- rivers:

  SpatRaster. Binary (0/1) raster for rivers.

- tri:

  SpatRaster. Continuous raster of topographic roughness (TRI). Used to
  increase cost in mountainous terrain.

- habitat:

  SpatRaster. Continuous raster of habitat suitability. Low values
  increase cost.

- w_ocean:

  Numeric. Weight applied to oceans (default 2000).

- w_lakes:

  Numeric. Weight applied to lakes (default 200).

- w_rivers:

  Numeric. Weight applied to rivers (default 20).

- w_tri:

  Numeric. Weight applied to TRI (default 1).

- w_habitat:

  Numeric. Weight applied to habitat suitability (default 1).

## Examples

``` r
if (FALSE) { # \dontrun{
# Prepare resistance raster

# this also can run internally in `population resistance`, 
# but for time sakes is best to prep ahead of time
# especially if treating multiple species in the same domain. 
res <- buildResistanceSurface(
  base_raster = base_rast,
  oceans = ocean_r,
  lakes = lakes_r,
  rivers = rivers_r,
  tri = tri_r
)
} # }
```
