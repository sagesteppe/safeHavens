## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)

## ----Install Package, message=FALSE, eval = F---------------------------------
# # remotes::install_github('sagesteppe/safeHavens')

## ----Load Libraries, message=FALSE,warnings=FALSE-----------------------------
library(safeHavens)
library(sf) ## vector operations
library(terra) ## raster operations
library(spData) ## basemap data
library(dplyr) ## general data handling
library(ggplot2) ## plotting 
library(patchwork) ## multiplots
set.seed(23) 

planar_proj <- "+proj=laea +lat_0=-15 +lon_0=-60 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

## ----Just Buffer it!----------------------------------------------------------
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

## ----echo = FALSE-------------------------------------------------------------
rm(x)

## ----Map themes---------------------------------------------------------------
x_extra_buff <- sf::st_buffer(x_buff, 100000) |> # add a buffer to 'frame' the maps
  sf::st_transform(4326)

americas <- spData::world
americas <- sf::st_crop(americas, sf::st_bbox(x_extra_buff)) |>
  dplyr::select(name_long)

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

## ----Point Based Sample, warning=FALSE, message = FALSE-----------------------
pbs <- PointBasedSample(x_buff, reps = 50, BS.reps = 333)
pbs.sf <- pbs[['Geometry']]

pbs.p <- map + 
  geom_sf(data = pbs.sf, aes(fill = factor(ID))) + 
#  geom_sf_label(data = pbs.sf, aes(label = ID), alpha = 0.4) + 
  labs(title = 'Point') + 
  coord_sf(expand = F)

pbs.p

## ----Equal Area Sample, message = FALSE, warning = FALSE----------------------
eas <- EqualAreaSample(x_buff, planar_proj = planar_proj) 

eas.p <- map + 
  geom_sf(data = eas[['Geometry']], aes(fill = factor(ID))) + 
#  geom_sf_label(data = eas.sf, aes(label = ID), alpha = 0.4) + 
  labs(title = 'Equal Area') + 
  coord_sf(expand = F)
eas.p

## ----echo = F-----------------------------------------------------------------
rm(eas, pbs, pbs.sf)

## ----Opportunistic Sample, message = FALSE------------------------------------
exist_pts <- sf::st_sample(x_buff, size = 10) |> 
   sf::st_as_sf() |> # ^^ randomly sampling 10 points in the species range
   dplyr::rename(geometry = x)

os <- OpportunisticSample(polygon = x_buff, n = 20, collections = exist_pts, reps = 50, BS.reps = 333)

os.p <- map + 
  geom_sf(data = os[['Geometry']], aes(fill = factor(ID))) + 
#  geom_sf_label(data = os.sf, aes(label = ID), alpha = 0.4) + 
  geom_sf(data = exist_pts, alpha = 0.4) + 
  labs(title = 'Opportunistic') + 
  coord_sf(expand = F)

os.p

## ----Isolatation by Distance Example, warning=FALSE, message = FALSE----------
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
ibdbs <- IBDBasedSample(
    x = v, 
    n = 20, 
    fixedClusters = TRUE, 
    template = predictors, 
    planar_proj = planar_proj
    )

ibdbs.p <- map + 
  geom_sf(data = ibdbs[['Geometry']], aes(fill = factor(ID))) + 
#  geom_sf_label(data = os.sf, aes(label = ID), alpha = 0.4) + 
  labs(title = 'IBDistance') + 
  coord_sf(expand = F)

## for the sake of comparing areas below, we will intersect this to the same extents as the earlier surfaces. 
ibdbs_crop <- sf::st_intersection(ibdbs[['Geometry']], sf::st_union(x_buff.sf))
ibdbs.p2 <- map + 
  geom_sf(data = ibdbs_crop, aes(fill = factor(ID))) + 
#  geom_sf_label(data = os.sf, aes(label = ID), alpha = 0.4) + 
  labs(title = 'IBDistance') + 
  coord_sf(expand = F)

ibdbs.p

## ----echo=F-------------------------------------------------------------------
rm(predictors, files, v, exist_pts, os)

## -----------------------------------------------------------------------------
ibr <- sf::st_read(
  file.path(system.file(package="safeHavens"), 'extdata', 'IBR.gpkg'), 
  quiet = TRUE)

ibr.p <- map + 
  geom_sf(data = ibr, aes(fill = factor(ID))) +
  labs(title = 'IBResistance') + 
  coord_sf(expand = F)


## ----Ecoregion Based Sample, message = F, warning = FALSE---------------------
neo_eco <- sf::st_read(
  file.path(system.file(package="safeHavens"), 'extdata', 'NeoTropicsEcoregions.gpkg'), 
  quiet = TRUE) |>
  dplyr::rename(geometry = geom)

head(st_drop_geometry(neo_eco)[,c('Provincias', 'Dominio', 'Subregion')]) |>
  knitr::kable()

x_buff <- sf::st_transform(x_buff, sf::st_crs(neo_eco))
ebs.sf <- PolygonBasedSample(x_buff, zones = neo_eco, n = 20, zone_key = 'Provincias')

# crop it to the other objects for plotting
ebs.sf <- st_crop(ebs.sf, bb)

ebs.p <- map + 
  geom_sf(data = ebs.sf , aes(fill = factor(allocation))) + 
  labs(title = 'Ecoregion') + 
  coord_sf(expand = F)

## ----echo = FALSE-------------------------------------------------------------
rm(bb, neo_eco, x_buff)

## ----Load data for Environmental Based Sample, message = F, results = 'hide'----
sdModel <- readRDS(
  file.path(system.file(package="safeHavens"), 'extdata',  'sdModel.rds')
  )

sdModel$RasterPredictions <- terra::unwrap(sdModel$RasterPredictions)
sdModel$Predictors <- terra::unwrap(sdModel$Predictors)
sdModel$PCNM <- terra::unwrap(sdModel$PCNM)

## -----------------------------------------------------------------------------
sdm <- terra::rast(
  file.path(system.file(package="safeHavens"), 'extdata',  'SDM_thresholds.tif')
  )
terra::plot(sdm)

## ----Environmental Based Sample-----------------------------------------------
rr <- RescaleRasters( # you may have already done this!
  model = sdModel$Model,
  predictors = sdModel$Predictors, 
  training_data = sdModel$TrainData, 
  pred_mat = sdModel$PredictMatrix
  )

# create a directory to hold the results from EBS real quick. 
# we will default to placing it in your current working directory. 
getwd() # this is where the folder is going to be created if you do not run the code below. 
p <- file.path(path.expand('~'), 'Documents') # in my case I'll dump it in Documents real quick, this should work on 

# optional, intentionally create a directory to hold results
# dir.create(file.path(p, 'safeHavens-Vignette'))

planar_proj <- "+proj=laea +lat_0=-15 +lon_0=-60 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

ENVIbs <- EnvironmentalBasedSample(
  pred_rescale = rr$RescaledPredictors, 
  write2disk = FALSE, # we are not writing, but showing how to provide some arguments
  path = file.path(p, 'safeHavens-Vignette'), 
  taxon = 'Bradypus_test', 
  f_rasts = sdm, 
  coord_wt = 2,
  n = 20, 
  lyr = 'Supplemented',
  fixedClusters = TRUE, 
  n_pts = 500, 
  planar_proj = planar_proj,
  buffer_d = 3,
  prop_split = 0.8
  )

## for the sake of comparing areas below, we will intersect this to the same extents as the earlier surfaces. 
ENVIbs_crop <- sf::st_intersection(ENVIbs[['Geometry']], sf::st_union(x_buff.sf))

ENVIbs.p <- map + 
  geom_sf(data = ENVIbs_crop, aes(fill = factor(ID))) + 
  #geom_sf_label(data = ENVIbs, aes(label = ID), alpha = 0.4) + 
  labs(title = 'Environmental') + 
  coord_sf(expand = FALSE)

ENVIbs.p

## ----Compare Environmental Based Samples, echo = F, message = F, warning = F----
ENVIbs_thresh <- EnvironmentalBasedSample(
  pred_rescale = rr$RescaledPredictors, 
  f_rasts = sdm, 
  planar_proj = planar_proj,
  lyr = 'Threshold'
  )

ENVIbs_clipped <- EnvironmentalBasedSample(
  pred_rescale = rr$RescaledPredictors, 
  f_rasts = sdm,
  planar_proj = planar_proj,
  lyr = 'Clipped'
  )

ENVIbs.p.thresh <- map + 
  geom_sf(data = ENVIbs_thresh[['Geometry']], aes(fill = factor(ID))) + 
  labs(title = 'Threshold') + 
  coord_sf(expand = FALSE)

ENVIbs.p.clipped <- map + 
  geom_sf(data = ENVIbs_clipped[['Geometry']], aes(fill = factor(ID))) + 
  labs(title = 'Clipped') + 
  coord_sf(expand = FALSE)

ENVIbs.p.supp <- map + 
  geom_sf(data = ENVIbs[['Geometry']], aes(fill = factor(ID))) + 
  labs(title = 'Supplemented') + 
  coord_sf(expand = FALSE)

ENVIbs.p.thresh + ENVIbs.p.clipped + ENVIbs.p.supp + 
  plot_layout(ncol = 3)

## ----echo = FALSE-------------------------------------------------------------
rm(ENVIbs_thresh, ENVIbs_clipped, ENVIbs.p.thresh, ENVIbs.p.clipped, p, ENVIbs, rr, ENVIbs.p.supp, sdm, map, sdModel)

## ----Plot all Sampling Schemes together---------------------------------------
pbs.p + eas.p + os.p  +  ibdbs.p2 + ibr.p + ENVIbs.p + 
  plot_layout(ncol = 3)

## ----echo = F-----------------------------------------------------------------
rm(pbs.p, eas.p, os.p, ibdbs.p, ebs.p, ENVIbs.p, planar_proj)

