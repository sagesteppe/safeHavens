## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)

## ----libraries, warning =F, message = F---------------------------------------
library(safeHavens)
library(sf) ## vector operations
library(terra) ## raster operations
library(dplyr) ## general data handling
library(ggplot2) ## plotting 

## ----import sample data-------------------------------------------------------
x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
x <- x[,c('lon', 'lat')]
x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)

planar_proj <- 3857 # Web Mercator for planar distance calcs

## ----clean vector data sets---------------------------------------------------
tri <- terra::rast(file.path(system.file(package = "safeHavens"), "extdata", "tri.tif"))
names(tri) <- 'tri'

rescale_rast <- function(r, new_min = 0, new_max = 1) {
  r_min <- global(r, "min", na.rm = TRUE)[[1]]
  r_max <- global(r, "max", na.rm = TRUE)[[1]]
  ((r - r_min) / (r_max - r_min)) * (new_max - new_min) + new_min
}

tri <- rescale_rast(tri, 0, 100)

lakes_v <- sf::st_read(
  file.path(system.file(package = "safeHavens"),  "extdata", "lakes.gpkg"),
    quiet = T) |>
  filter(TYPE == 'Lake') |> # remove reservoirs for our purposes (time scale)
  select(geometry = geom) |> 
  mutate(Value = 1) |> # this will be rasterized. 
  terra::vect()

ocean_v <- sf::st_read(
  file.path(system.file(package = "safeHavens"), "extdata", "oceans.gpkg"),
   quiet = T) |>
  select(geometry = geom) |>
  mutate(Value = 1) |>
  terra::vect()

rivers_v <- sf::st_read(
  file.path(system.file(package = "safeHavens"), "extdata", "rivers.gpkg"),
   quiet = T) |>
  select(geometry = geom) |>
  mutate(Value = 1) |>
  terra::vect()


## ----crop data and project----------------------------------------------------
x_buff <- sf::st_transform(x, planar_proj) |>
  # huge buffer for the bbox. 
  st_buffer(200000) |> 
  st_transform(crs(lakes_v)) |>
  st_as_sfc() |> 
  st_union() |>
  vect() |>
  ext()

lakes_v <- crop(lakes_v, x_buff)
rivers_v <- crop(rivers_v, x_buff)
ocean_v <- crop(ocean_v, x_buff)
tri <- crop(tri, x_buff)

## ----rasterize data, message=F,warning=F--------------------------------------
ocean_r <- rasterize(ocean_v, tri, field = 'Value', background = 0.1)
lakes_r <- rasterize(lakes_v, tri, field = 'Value',  background = 0.1)
rivers_r <- rasterize(rivers_v, tri, field = 'Value',  background = 0.1)

par(mfrow=c(2, 2))
plot(rivers_r)
plot(lakes_r)
plot(ocean_r)
plot(tri)
par(mfrow=c(1,1))

rm(ocean_v, lakes_v, rivers_v)

## ----create resistance surface------------------------------------------------
res_surface <- buildResistanceSurface(
  base_raster = rast(ocean_r), # template

  oceans = ocean_r, # rasters
  lakes = lakes_r,
  rivers = rivers_r,
  tri = tri,

  w_ocean = 120, # weights to multiple  input raster values by -- 
  w_lakes = 50, # is 1 * 50
  w_rivers = 20, # is 1 *20
  w_tri = 4 # ranges from 1~30, so from (1-30)*4 up to a weight of ~120. 
)

plot(res_surface)

## ----generate population graphs and distance matrix, warning=F----------------
pop_res_graphs <- populationResistance(
  base_raster = rast(ocean_r),
  populations_sf = x,
  n = 150,
  planar_proj = 3857,
  buffer_dist = 125000,
  resistance_surface = res_surface,
  graph_method = 'complete'
)

names(pop_res_graphs)

## ----plot the buffered populations--------------------------------------------
plot(pop_res_graphs$pop_raster)

## ----run IBRBasedSample, warning = F------------------------------------------
ibr_groups <- IBRSurface(
  base_raster = rast(ocean_r),
  resistance_surface = res_surface, 
  pop_raster = pop_res_graphs$pop_raster,
  pts_sf = pop_res_graphs$sampled_points, 
  ibr_matrix = pop_res_graphs$ibr_matrix, 
  fixedClusters = TRUE,
  n = 20,
  planar_proj = 3857
)


## ----Plot cluster classified points-------------------------------------------
classified_pts = ibr_groups$points
plot(res_surface)
points(
  classified_pts, 
  pch = as.numeric(as.factor(classified_pts$ID)),  # different symbol per ID
  col = rainbow(length(unique(classified_pts$ID)))[as.factor(classified_pts$ID)],
  cex = 2,
  lwd = 3
  )

## ----plot the classified clusters---------------------------------------------
ggplot() +
  geom_sf(data = ibr_groups$geometry, aes(fill = factor(ID))) + 
  theme_minimal()

## ----warning = F--------------------------------------------------------------
out <- PolygonBasedSample(
  x = st_union(ibr_groups$geometry),
  zones = ibr_groups$geometry, 
  zone_key = 'ID', 
  n  = 30
  )

ggplot(data = out) + 
  geom_sf(aes(fill = factor(allocation))) + 
  theme_minimal() + 
  labs(fill = 'Ct. collections\nper cluster') + 
  theme(legend.position = 'bottom')

## ----save results for getting started final plot, eval = F, echo = F----------
# #p <- file.path('~', 'Documents', 'assoRted', 'safeHavens', 'inst', 'extdata')
# #sf::st_write(ibr_groups$geometry, file.path(p, 'IBR.gpkg'), driver = 'gpkg')

## -----------------------------------------------------------------------------
sdm <- terra::rast(
  file.path(system.file(package="safeHavens"), 'extdata',  'SDM_thresholds.tif')
  )
plot(sdm['Predictions'])

## -----------------------------------------------------------------------------
inverted_sdm <- 1 - sdm['Predictions']
plot(inverted_sdm)

## -----------------------------------------------------------------------------
inverted_sdm <- terra::ifel(inverted_sdm < 0.5, 0.01, inverted_sdm)
plot(inverted_sdm)

## -----------------------------------------------------------------------------
inverted_sdm = round(inverted_sdm * 100, 0)
plot(inverted_sdm)

## -----------------------------------------------------------------------------
inverted_sdm <- crop(inverted_sdm, ocean_r)
inverted_sdm <- resample(inverted_sdm, ocean_r)
inverted_sdm[is.na(inverted_sdm)] <- 99

plot(inverted_sdm)

## ----plot non-linear transformations------------------------------------------
x <- 1:100
rs <- function(x, to = c(0, 100)) {
  rng <- range(x, na.rm = TRUE)
  to[1] + (x - rng[1]) / diff(rng) * diff(to)
}

# display some transformations. 
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
plot(x, rs(sqrt(x)), type = "l", main = "sqrt(x)")
plot(x, rs(log1p(x)), type = "l", main = "log1p(x)")
plot(x, rs(asinh(x)), type = "l", main = "asinh(x)")
plot(x, rs(x^2), type = "l", main = "x^2: polynomial")
plot(x, rs(x^3), type = "l", main = "x^3: polynomial")
plot(x, rs(x^4), type = "l", main = "x^4 polynomial")

## ----remove var,echo=F--------------------------------------------------------
rm(x)

## ----rescale TRI--------------------------------------------------------------
tri_rscl <- app(tri, function(x) { rs(x^2, to = c(1, 80)) })

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(tri, main = 'Original')
plot(tri_rscl, main = 'Rescaled')

## ----build a polynomial surface-----------------------------------------------
res_surface <- buildResistanceSurface(

  base_raster = rast(ocean_r), # template

  oceans = ocean_r, # rasters
  lakes = lakes_r,
  rivers = rivers_r,
  #tri = tri, # dont use me! i rescale. 
  habitat = inverted_sdm,
  addtl_r = tri_rscl,

  w_ocean = 120, # weights to multiple  input raster values by -- 
  w_lakes = 70, # is 1 * 50
  w_rivers = 40, # is 1 *20
  w_habitat = 0.5,
  #w_tri = 4 # don't use me! use addtl_wt
  addtl_w = 1
)
par(mfrow=c(1,1))
plot(res_surface)

