# Create hexagonal grid based polygons over a geographic range

This function creates 20 grid cells over a geographic area (`x`),
typically a species range.

## Usage

``` r
GridBasedSample(x, planar_proj, gridDimensions)
```

## Arguments

- x:

  An SF object or terra spatraster. the range over which to generate the
  clusters.

- planar_proj:

  Numeric, or character vector. An EPSG code, or a proj4 string, for a
  planar coordinate projection, in meters, for use with the function.
  For species with very narrow ranges a UTM zone may be best (e.g. 32611
  for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). Otherwise
  a continental scale projection like 5070 See
  https://projectionwizard.org/ for more information on CRS. The value
  is simply passed to sf::st_transform if you need to experiment.

- gridDimensions:

  A single row form the ouput of `TestGridSizes` with the optimal number
  of grids to generate.

## Value

An simple features (sf) object containing the final grids for saving to
computer. See the vignette for questions about saving the two main types
of spatial data models (vector - used here, and raster).

## Examples

``` r
if (FALSE)  # not ran to bypass CRAN check time limits. ~6 seconds to treat Rhode Island. 
ri <- spData::us_states |> 
dplyr::filter(NAME == 'Rhode Island') |>
  sf::st_transform(32615)

sizeOptions <- TestGridSizes(ri)
#> Error: object 'ri' not found
head(sizeOptions) # in this case let's shoot for 33 and see what happens
#> Error: object 'sizeOptions' not found
sizeOptions <- sizeOptions[sizeOptions$Name == 'Original',]
#> Error: object 'sizeOptions' not found

output <- GridBasedSample(ri, 5070, gridDimensions = sizeOptions)
#> Error: object 'ri' not found
plot(output$Geometry)
#> Error: object 'output' not found
 # \dontrun{}
```
