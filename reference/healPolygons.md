# Clean up unioned geometries - part 1

this function uses sf::st_snap to remove small lines and other artifacts
associated with the unioning of polygons. This is ran within `snapGrids`

## Usage

``` r
healPolygons(x)
```

## Arguments

- x:

  most of the output of `snapgrids`
