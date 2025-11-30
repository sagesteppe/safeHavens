# Getting Started

### set up

`safeHavens` can be installed directly from github.

``` r
remotes::install_github('sagesteppe/safeHavens') 
```

``` r
library(safeHavens)
library(ggplot2)
library(sf)
library(terra)
library(spData)
library(dplyr)
library(patchwork)
set.seed(23) 
planar_proj <- '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs'
```

### Defining a Species Range or Domain for Sampling

Central to the sampling schemes in `safeHavens` is a species range or
domain for sampling. For example, depending on the goals of the
collection, a curator may want to sample across the entire range of a
species. Alternatively one may be interested in sampling only a portion
of the range, e.g. a country, state, or ecoregion. Either of these
scenarios can be accomplished with the package. Here we will show how to
create a species range from occurrence data, and then use that range to
run the various sampling schemes.

Below we will use `sf` to simply buffer occurrence points to create a
species range across multiple South American nations.

``` r
x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
x <- x[,c('lon', 'lat')]
x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)
x_buff <- sf::st_transform(x, planar_proj) |>
  # we are working in planar metric coordinates, we are
  # buffer by this many / 1000 kilometers. 
  sf::st_buffer(125000) |> 
  sf::st_as_sfc() |> 
  sf::st_union()

plot(x_buff)
```

![](GettingStarted_files/figure-html/Just%20Buffer%20it!-1.png)

Alternatives to this include a simple convex hull around a species, if
it is widespread throughout an area, or masking a binary SDM surface as
the domain.

## Prep a map background

We will use the `spData` package which uses
[naturalearth](https://www.naturalearthdata.com/) data for it’s `world`
data and is suitable for creating effective maps at a variety of
resolutions.

``` r
x_extra_buff <- sf::st_buffer(x_buff, 100000) |> # add a buffer to 'frame' the maps
  sf::st_transform(4326)

americas <- spData::world
americas <- sf::st_crop(americas, sf::st_bbox(x_extra_buff)) |>
  dplyr::select(name_long)
#> Warning: attribute variables are assumed to be spatially constant throughout
#> all geometries

bb <- sf::st_bbox(x_extra_buff)

map <- ggplot() + 
  geom_sf(data = americas) + 
  theme(
    legend.position = 'none', 
    panel.background = element_rect(fill = "aliceblue"), 
    panel.grid.minor.x = element_line(colour = "red", linetype = 3, linewidth  = 0.5), 
    axis.ticks=element_blank(),
    axis.text=element_blank(),
    plot.background=element_rect(colour="steelblue"),
    plot.margin=grid::unit(c(0,0,0,0),"cm"),
    axis.ticks.length = unit(0, "pt"))+ 
  coord_sf(xlim = c(bb[1], bb[3]), ylim = c(bb[2], bb[4]), expand = FALSE)

rm(x_extra_buff, americas)
```

## Running the Various Sample Design Algorithms

Now that we have some data which can represent species ranges, we can
run the various sampling approaches. The table in the introduction is
reproduced here.

| Function                   | Description                             | Comp. | Envi. |
|----------------------------|-----------------------------------------|-------|-------|
| `GridBasedSample`          | Creates and merges *n* grids over area  | L     | L     |
| `PointBasedSample`         | Creates points to make pieces over area | L     | L     |
| `EqualAreaSample`          | Breaks area into similar size pieces    | L     | L     |
| `OpportunisticSample`      | Using PBS with existing records         | L     | L     |
| `IBDBasedSample`           | Breaks species range into clusters      | H     | M     |
| `EcoregionBasedSample`     | Using existing ecoregions to sample     | L     | H     |
| `EnvironmentalBasedSample` | Uses correlations from SDM to sample    | H     | H     |

*Note in this table ‘Comp.’ and ‘Envi.’ refer to computational and
environmental complexity respectively, and range from low (L) through
medium to high.*

### Grid Based Sample

Grids are very useful for sampling contiguous things. Species ranges are
often not contiguous; however curators and analysts in our
geographically grand ecosystems, e.g. the steppes, prairies, tundra, and
taiga might find these useful.

You can look at the output below to see why grids are not great for
*this* type of problem.

The first step in grid sampling is determining an OK number of grids to
try and draw as a starting point, if we want 20 collections we will need
more than 20 grids, because several will be merged into the larger ones.
Using the aspect ratio of a simple bounding box around the area will be
analyzing, the `TestGridSizes` function will determine a default number
of grids (‘Original’) for testing. Using these defaults it will create a
few other sets of grids as well, by either removing one of two grids per
direction. Theoretically you could automate grid selection by comparing
the number of grids and the minimization of variance. To be safe I
wouldn’t consider configurations which generate less than 25 of these
initial grids.

``` r
tgs <- TestGridSizes(x_buff)
print(tgs)
#>       Name Grids  Variance GridNOx GridNOy
#> 1 Smallest    34  599.9513       8       6
#> 2  Smaller    28  813.2043       7       5
#> 3 Original    24  862.5458       6       4
#> 4   Larger    18  681.7667       5       3
#> 5  Largest    13 1182.8457       4       2

plot(tgs$Grids, tgs$Variance, xlab = 'Grid Number', ylab = 'Variance',
     main = 'Number of grids and areas overlapping species range')
text(tgs$Grids, tgs$Variance + 25, labels=tgs$Name, cex= 0.7)
abline(v=20, col="red")
abline(v=25, col="orange")
```

![](GettingStarted_files/figure-html/Test%20Grid%20Sizes-1.png)

``` r

tgs <- tgs[tgs$Name=='Smaller',]
```

Essentially we need more than 20 grids, but realistically (albeit from
limited informal testing) using more than 25 grids - depending on the
complexity of the species range - tends to be an effective floor. In the
table and plot above I opt for using the ‘Smaller’ option, with 28 grids
generated by prompting
[`sf::st_make_grid`](https://r-spatial.github.io/sf/reference/st_make_grid.html)
with 7 grids in the X direction and 5 in the Y direction. You can *kind
of* think about this like an elbow plot, but the samples are so few you
won’t get the characteristic shape.

``` r
grid_buff <- GridBasedSample(x_buff, planar_proj, gridDimensions = tgs) 

gbs.p <- map + 
  geom_sf(data = grid_buff, aes(fill = factor(ID))) + 
 # geom_sf_label(data = grid_buff, aes(label = Assigned), alpha = 0.4) +  # on your computer, doesnt work at vignette size
  labs(title = 'Grids')  + 
  coord_sf(expand = F)

gbs.p
```

![](GettingStarted_files/figure-html/Grid%20Based%20Samples-1.png)

### Point Based Sample

With the grids we drew the pre-specified number of grids across the
species range, and then merged them together as required to get the
results. We will essentially do the inverse in this step, rather than
drawing boundaries - i.e. grid cells, we will draw the centers. This
essentially allows the features to ‘grow’ a little more naturally. I
also think these results work a little bit better on a fragmented range,
their is still some odd clipping, where minor portions of a section of
range are assigned to a different grid, but in general a little bit
better.

``` r
pbs <- PointBasedSample(x_buff)
pbs.sf <- pbs$Geometry

st_is_valid(x_buff)
#> [1] TRUE

pbs.p <- map + 
  geom_sf(data = pbs.sf, aes(fill = factor(ID))) + 
#  geom_sf_label(data = pbs.sf, aes(label = ID), alpha = 0.4) + 
  labs(title = 'Point') + 
  coord_sf(expand = F)
pbs.p
```

![](GettingStarted_files/figure-html/Point%20Based%20Sample-1.png)

### Equal Area Sample

Perhaps the simplest method which is offered in `safeHavens` is
`EqualAreaSample`. It simply creates many points, `pts` defaulting to
5000, within our target domain and subjects them to k-means sampling
where the groups are specified by `n`, our target number of collections.
The individual points assigned to each group are merged into polygons
which ‘take’ up all of the geographic space, and are intersected back to
the species range, and the area of each polygon is then measured. This
process will be ran a few times, defaulting to 100 `reps`, and the set
of polygons which was created during these reps with the smallest
variance in polygon size will be selected and returned.

This differs from point based sampling in that the above instance, we
start with a few regularly spaced points to grow from, here we take a
step back and by using many points let the clusters grow themselves to
similar sizes.

``` r
eas <- EqualAreaSample(x_buff, planar_projection = planar_proj) 
#> Warning: did not converge in 10 iterations

eas.p <- map + 
  geom_sf(data = eas$Geometry, aes(fill = factor(ID))) + 
#  geom_sf_label(data = eas.sf, aes(label = ID), alpha = 0.4) + 
  labs(title = 'Equal Area') + 
  coord_sf(expand = F)
eas.p
```

![](GettingStarted_files/figure-html/Equal%20Area%20Sample-1.png)

The results look quite similar to point based sample.

### Opportunistic Sample

Many curators are interested in how much they can embed their existing
collections into a sampling framework. The function
`OpportunisticSample` makes a few minor modifications to the point based
sample to maximize an existing collection. It doesn’t always work
exceptionally, especially when a couple collections are very close to
each other, but it may be a beneficial tool to have in the belt. As we
have observed, the three previous sampling schemes end of with somewhat
similar results - so we took used the `PointBasedSample` as the
framework we embedded this function into it.

Essentially this combines the approach of point based sampling, but
forces that the clusters are based around the existing accessions. It
attempts to ‘center’ the existing collections within clusters, but this
can be nearly impossible for a variety of reasons.

``` r
exist_pts <- sf::st_sample(x_buff, size = 10) |> 
   sf::st_as_sf() |> # ^^ just randomly sampling 10 points in the species range
   dplyr::rename(geometry = x)

os <- OpportunisticSample(polygon = x_buff, n = 20, collections = exist_pts)

os.p <- map + 
  geom_sf(data = os$Geometry, aes(fill = factor(ID))) + 
#  geom_sf_label(data = os.sf, aes(label = ID), alpha = 0.4) + 
  geom_sf(data = exist_pts, alpha = 0.4) + 
  labs(title = 'Opportunistic') + 
  coord_sf(expand = F)

os.p
```

![](GettingStarted_files/figure-html/Opportunistic%20Sample-1.png)

Here, the grids have been aligned around the points. This can lead to
some oddly shaped clusters, but *a bird in hand is worth two in the
bush.*

### Isolation by Distance Based Sample

Isolation by Distance is the fundamental idea behind this package. This
function explicitly uses IBD to develop a sampling scheme, and does not
obfuscate it with any other parameters.

Note that this function requires a raster input, rather than a vector.

``` r
files <- list.files( 
  path = file.path(system.file(package="dismo"), 'ex'),
  pattern = 'grd',  full.names=TRUE ) 
predictors <- terra::rast(files) 

x_buff.sf <- sf::st_as_sf(x_buff) |> 
  dplyr::mutate(Range = 1) |> 
  sf::st_transform( terra::crs(predictors))

# and here we specify the field/column with our variable we want to become an attribute of our raster
v <- terra::rasterize(x_buff.sf, predictors, field = 'Range') 

# now we run the function demanding 20 areas to make accessions from, 
ibdbs <- IBDBasedSample(x = v, n = 20, fixedClusters = TRUE, template = predictors)

ibdbs.p <- map + 
  geom_sf(data = ibdbs, aes(fill = factor(ID))) + 
#  geom_sf_label(data = os.sf, aes(label = ID), alpha = 0.4) + 
  labs(title = 'IBD') + 
  coord_sf(expand = F)
ibdbs.p
```

![](GettingStarted_files/figure-html/Isolatation%20by%20Distance%20Example-1.png)

``` r

rm(predictors, files, v, x_buff.sf, exist_pts, os)
```

Because these data were processed from a raster, they have lienar edges,
representing raster tiles. However, it is evident that the borders of
the clusters are more natural looking than in the previous (and future)
sampling schemes.

### Ecoregion Based Sample

This is the most commonly implemented method for guiding native seed
collection in North America. However, I am not sure exactly how
practitioners all implement it, and whether the formats of application
are consistent among practitioners! For these reasons a few different
sets of options are supported for a user.

For general usage, two parameters are always required `x` which is the
species range as an sf object, and `ecoregions`, the sf object
containing the ecoregions of interest. The `ecoregions` file does not
need to be subset to the range of `x` quite yet - the function will take
care of that. Additional arguments to the function include as usual `n`
to specify how many accession we are looking for in our collection. Two
additional arguments relate to whether we are using Omernik Level 4
ecoregions data or ecoregions (or biogeographic regions) from another
source. These are `OmernikEPA`, and `ecoregion_col`, if you are using
the official EPA release of ecoregions then both of these are optional,
however if you are not using the EPA product than both *should* be
supplied - but only the `ecoregion_col` argument is totally necessary.
This column should contain unique names for the highest resolution level
ecoregion you want to use from the data set, for many data sets, such as
our example we call ‘neo_eco’ this may be the only field with ecolevel
information!

``` r
neo_eco <- sf::st_read(
  file.path(system.file(package="safeHavens"), 'extdata', 'NeoTropicsEcoregions.gpkg'), 
  quiet = TRUE) |>
  dplyr::rename(geometry = geom)
head(neo_eco[,c(1, 3, 4, 6, 11)])
#> Simple feature collection with 6 features and 4 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -103.0432 ymin: -31.25308 xmax: -34.79344 ymax: 26.91751
#> Geodetic CRS:  WGS 84
#>                  Provincias      Region      Dominio
#> 1 Araucaria Forest province Neotropical       Parana
#> 2          Atacama province Neotropical         <NA>
#> 3         Atlantic province Neotropical       Parana
#> 4           Bahama province Neotropical         <NA>
#> 5     Balsas Basin province Neotropical Mesoamerican
#> 6         Caatinga province Neotropical      Chacoan
#>                        Subregion                       geometry
#> 1                        Chacoan MULTIPOLYGON (((-53.58012 -...
#> 2 South American Transition Zone MULTIPOLYGON (((-69.42981 -...
#> 3                        Chacoan MULTIPOLYGON (((-48.41217 -...
#> 4                      Antillean MULTIPOLYGON (((-77.58593 2...
#> 5                      Brazilian MULTIPOLYGON (((-97.37265 1...
#> 6                        Chacoan MULTIPOLYGON (((-35.56652 -...

x_buff <- sf::st_transform(x_buff, sf::st_crs(neo_eco))
ebs.sf <- EcoregionBasedSample(x_buff, neo_eco, OmernikEPA = FALSE, ecoregion_col = 'Provincias')
#> Warning: attribute variables are assumed to be spatially constant throughout
#> all geometries
#> Warning: attribute variables are assumed to be spatially constant throughout
#> all geometries

# for plotting let's crop it to the other objects
ebs.sf <- st_crop(ebs.sf, bb)
#> Warning: attribute variables are assumed to be spatially constant throughout
#> all geometries

ebs.p <- map + 
  geom_sf(data = ebs.sf, aes(fill = factor(n))) + 
  labs(title = 'Ecoregion') + 
  coord_sf(expand = F)
ebs.p
```

![](GettingStarted_files/figure-html/Ecoregion%20Based%20Sample-1.png)

This output differs from the others we will see, here we have depicted
the number of collections to be made per ecoregion. Because the number
of ecoregions is greater than our requested sample size, the return
object can only take on two values - no collections, or one collection.

### Environmental Based Sample

The environmental based sample can only be conducted if you have the
species distribution model data. Included in the data directory of the
folder are all of the objects required to run this example for the
species. We will load them here.

## load the SDM predictions

The package also has an SDM prediction saved in the data we can just
load that too for a couple comparisons.

``` r
sdm <- terra::rast(file.path(system.file(package="safeHavens"),  'extdata', 'Bradypus_test.tif'))
terra::plot(sdm)
```

![](GettingStarted_files/figure-html/unnamed-chunk-5-1.png)

``` r
sdModel <- readRDS(
  file.path(system.file(package="safeHavens"), 'extdata',  'sdModel.rds')
  )

sdModel$Predictors <- terra::rast(
  file.path(system.file(package="safeHavens"), 'extdata', 'Predictors.tif')
)
```

Once these data are loaded into R, we will scale the rasters (using
`RescaleRasters`) which will serve as surfaces to predict from (this is
also done above!), then we will run the algorithm
(`EnvironmentalBasedSample`). However, before we run the algorithm we
will need to create a directory (also called a ‘folder’), on our
computers to save the results from the function
`EnvironmentalBasedSample`. Whereas earlier in this vignette we
showcased that the functions generated the species distribution model,
and us saving the results were a two stage process (e.g. to create the
SDM and associated products we used: `elasticSDM`, `PostProcessSDM`, and
`RescaleRasters`, before finally saving relevant data with
`writeSDMresults`), this function produces both the product and writes
out ancillary data simultaneously.  
This approach was chosen as this function is only writing out four
objects: 1) the groups as vector data, 2) and the groups as raster data,
3) the k-nearest neighbors (knn) model used to generate these clusters,
and 4) the confusion matrix associated with testing the knn model.

``` r
rr <- RescaleRasters( # you may have already done this!
  model = sdModel$Model,
  predictors = sdModel$Predictors, 
  training_data = sdModel$TrainData, 
  pred_mat = sdModel$PredictMatrix)

# create a directory to hold the results from EBS real quick. 
# we will default to placing it in your current working directory. 
# If you are a data management freak don't worry too much about this. 
# The code to remove the directory will be provided below. 
getwd() # this is where the folder is going to located IF YOU DON'T RUN the code below. 
#> [1] "/home/runner/work/safeHavens/safeHavens/vignettes"
p <- file.path('~', 'Documents') # in my case I'll dump it in Documents real quick, this should work on 
# Linux and Mac, but I don't think Windows? 
# dir.create(file.path(p, 'safeHavens-Vignette')) # now create the directory. 

ENVIbs <- EnvironmentalBasedSample(
  pred_rescale = rr$RescaledPredictors, 
  write2disk = FALSE, 
  path = file.path(p, 'safeHavens-Vignette'), # we are not writing, but showing how to provide argument
  taxon = 'Bradypus_test', 
  f_rasts = sdm, n = 20, 
  lyr = 'Supplemented',
  fixedClusters = TRUE, 
  n_pts = 500, 
  planar_projection = planar_proj,
  buffer_d = 3, prop_split = 0.8)
#> Joining with `by = join_by(x, y)`
#> Warning in st_point_on_surface.sfc(st_geometry(x)): st_point_on_surface may not
#> give correct results for longitude/latitude data

ENVIbs.p <- map + 
  geom_sf(data = ENVIbs, aes(fill = factor(ID))) + 
  #geom_sf_label(data = ENVIbs, aes(label = ID), alpha = 0.4) + 
  labs(title = 'Environmental') + 
  coord_sf(expand = FALSE)
#> Coordinate system already present.
#> ℹ Adding new coordinate system, which will replace the existing one.

ENVIbs.p
```

![](GettingStarted_files/figure-html/Environmental%20Based%20Sample-1.png)

The function `EnvironmentalBasedSample` can take any of the three binary
rasters created by `PostProcessSDM` as arguments for the template. Here
we showcase the different results from using each of them.

![](GettingStarted_files/figure-html/Compare%20Environmental%20Based%20Samples-1.png)

These plots are able to showcase the difference in results depending on
which of the three input rasters are utilized. As with all of the
sampling schemes, results vary widely based on the spatial extents which
the functions are applied to. Using the SDM output which have undergone
thresholding results in the largest classified area. At first glance the
results may seem very different, but if you look at central america,
they are largely consistent, as they are near the Andes; large
differences do exist in the Amazon Basin, but even there some alignment
between the systems is evident. Accordingly, the surface used for a
species should match some evaluation criterion.

Using the threshold raster surface is a very good option if we do not
want to ‘miss’ too many areas, whereas the clipped and supplemented
options may be better suited for scenarios where we do not want to draw
up clusters, which lack any populations which can be collected from.

## Comparision of different sampling schemes

So, we’ve made got some maps for you to look at! They all look
relatively similar to me when plotted one after another, let’s plot them
all simultaneously and see if that’s still the case.

``` r
gbs.p + pbs.p + eas.p + os.p  +  ibdbs.p + ebs.p + ENVIbs.p + 
  plot_layout(ncol = 3)
```

![](GettingStarted_files/figure-html/Plot%20all%20Sampling%20Schemes%20together-1.png)

Again, the top three figures appear quite similar, with the
`Opportunistic` method only deviating slightly form them. In my mind
isolation by distance (IBD) show the biggest different, it seems to have
made the most *sense* of the naturally occurring patchiness of the
species range. Ecoregion SEEMS…. Environmental also seems to partition
the feature space quite well. Notably drawing a couple clusters in the
Pacific lowlands and Northern Andes mountains.
