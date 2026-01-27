library(terra)
library(sf)
library(tidyverse)
library(safeHavens)

x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
x <- x[,c('lon', 'lat')]
x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)

planar_proj <- 3857 

setwd(file.path('~', 'Documents', 'assoRted', 'safeHavens-data'))

sf_data <- file.path('~', 'Documents', 'assoRted', 'safeHavens-data')
tri <- terra::rast(file.path(sf_data, 'tri_5KMmd_GMTEDmd.tif'))
names(tri) <- 'tri'



lakes_v <- sf::st_read(
  file.path(sf_data, 'GLWD-level1', 'glwd_1.shp'), quiet = T) |>
  st_set_crs(4326) |> ## messed up shape, manually tell it the crs. 
  select(TYPE, geometry) |> 
  filter(TYPE == 'Lake') |> # remove reservoirs for our purposes (time scale)
  mutate(Value = 1) |> # this will be rasterized. 
  terra::vect()

ocean_v <- sf::st_read(
  file.path(sf_data, 'ne_10m_ocean', 'ne_10m_ocean.shp'), quiet = T) |>
  select(geometry) |>
  mutate(Value = 1) |>
  terra::vect()

rivers_v <- sf::st_read(
  file.path(sf_data, 'rivers_world_47950', 'rivers_world_47950.shp'), quiet = T) |>
  filter(Strahler > 3) |> # just showing this option exists... 
  select(geometry) |>
  mutate(Value = 1) |>
  terra::vect()


x_buff <- sf::st_transform(x, planar_proj) |>
  # huge buffer for the bbox. 
  st_buffer(200000) |> 
  st_transform(crs(lakes_v)) |>
  st_as_sfc() |> 
  st_union() |>
  vect() |>
  ext()

lakes_v <- crop(lakes_v, x_buff)  |>
  st_as_sf() |>
  st_make_valid() |>
  st_simplify() |>
  st_make_valid()
rivers_v <- crop(rivers_v, x_buff)|>
  st_as_sf() |>
  st_simplify() |>
  st_make_valid()
ocean_v <- crop(ocean_v, x_buff)|>
  st_as_sf() |>
  st_simplify() |>
  st_make_valid()
tri <- crop(tri, x_buff)

p <- file.path('~', 'Documents', 'assoRted', 'safeHavens', 'inst', 'extdata')
sf::st_write(lakes_v, file.path(p, 'lakes.gpkg'), driver = 'gpkg')

sf::st_write(rivers_v, file.path(p, 'rivers.gpkg'), driver = 'gpkg')
sf::st_write(ocean_v, file.path(p, 'oceans.gpkg'), driver = 'gpkg')

tri <- as.int(tri)
terra::writeRaster(tri, file.path(p, 'tri.tif'))
