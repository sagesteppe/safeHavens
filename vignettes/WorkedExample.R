## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)

## ----load packages, warning = F, message = F----------------------------------
library(safeHavens)
library(sf) # spatial data
library(rgbif) # species occurrence data.
library(spData) # example cartography data
library(dplyr) # general data handling
library(tidyr) # general data handling
library(purrr) # sub for lapplying functions across lists
library(ggplot2) ## for maps 

## ----download species occurrence data from gbif-------------------------------
## small subset of useful columns for example
cols = c('decimalLatitude', 'decimalLongitude', 'dateIdentified', 'species', 'acceptedScientificName', 'datasetName', 
  'coordinateUncertaintyInMeters', 'basisOfRecord', 'institutionCode', 'catalogNumber')

## download species data using scientificName, can use keys and lookup tables for automating many taxa. 
cymu <- rgbif::occ_search(scientificName = "Vesper multinervatus", limit = 1000)

### check to see what CRS are in here, these days usually standardized to wgs84 (epsg:4326)
table( cymu[['data']]['geodeticDatum']) 

## subset the data to relevant columns 
cymu_cols <- cymu[['data']][,cols]

## repeat this again so a second set of data are on hand 
bowa <- rgbif::occ_search(scientificName = "Bouteloua warnockii", limit = 1000)
bowa_cols <- bowa[['data']][,cols]

## contrived multispecies example. 
spp <- bind_rows(bowa_cols, cymu_cols) |>
  drop_na(decimalLatitude, decimalLongitude) |> # any missing coords need dropped. 
  st_as_sf(coords = c( 'decimalLongitude', 'decimalLatitude'), crs = 4326, remove = F)

rm(cymu, bowa, cols, cymu_cols, bowa_cols)

## ----check gbif data----------------------------------------------------------
western_states <- spData::us_states |> ## for making a quick basemap. 
  dplyr::filter(REGION == 'West' & ! NAME %in% c('Montana', 'Washington', 'Idaho', 'Oregon', 'Wyoming') |
   NAME %in% c('Oklahoma', 'Texas', 'Kansas')) |>
  dplyr::select(NAME, geometry) |>
  st_transform(4326)

ggplot() +
  geom_sf(data = western_states) +
  geom_sf(data = spp, aes(color = species, shape = species)) +
  theme_void() +
  theme(legend.position = 'bottom')

## check the outlying record. 
arrange(spp, by = decimalLatitude, desc=FALSE) |>
  head(5) |>
  knitr::kable()

## remove it based on it's latitude. 
spp <- filter(spp, decimalLatitude <= 40)

## ----create basemap, warning=FALSE--------------------------------------------
bb <- st_transform(spp, 5070) |>
  st_buffer(100000) |>
  st_transform(4326) |>
  st_bbox()

western_states <- st_crop(western_states, bb)

base <- ggplot() +
  geom_sf(data = western_states, color = 'white') +
  geom_sf(data = spp, aes(color = species, fill = species)) +
  theme_void() +
  theme(legend.position = 'bottom')

base

## ----remove western states, echo = F------------------------------------------
rm(western_states)

## ----species ranges-----------------------------------------------------------
sppL <- split(spp, f = spp$species)

concavities <- function(x, d, rat){

  species <- x[['species']][1]
  out <- st_transform(x, 5070) |>
    st_buffer(dist = d) |>
    st_union() |> 
    st_concave_hull(ratio = rat) |> 
    st_sf() |>
    rename(geometry = 1) |>
    mutate(species, .before = geometry)

}

spp_concave <- sppL |>
  purrr::map(~ concavities(.x, d = 20000, rat = 0.4))
spp_convex <- sppL |>
  purrr::map(~ concavities(.x, d = 20000, rat = 1.0))

rm(concavities)

## ----plot concave hulls-------------------------------------------------------
base +
  geom_sf(data = spp_concave[[1]], aes(fill = species), alpha = 0.2) +
  geom_sf(data = spp_concave[[2]],  aes(fill = species), alpha = 0.2) 

## ----plot convex hulls--------------------------------------------------------
base +
  geom_sf(data = spp_convex[[1]],  aes(fill = species), alpha = 0.2) +
  geom_sf(data = spp_convex[[2]], aes(fill = species), alpha = 0.2) 

rm(spp_convex)

## ----determine radius for buffering points, eval = F, echo = F----------------
# buffered <- sppL |>
#   purrr::map(~ vario_estimate(.x, bb_dist = 50000))
# 
# lapply(buffered, `[[`, 'range')
# 
# #buffered +
# #  geom_sf(data = buffered[[1]],  aes(fill = species), alpha = 0.2) +
# #  geom_sf(data = buffered[[2]], aes(fill = species), alpha = 0.2)

## ----perform sampling---------------------------------------------------------
eas <- spp_concave |>
  purrr::map(~ EqualAreaSample(.x, n = 10, pts = 250, planar_proj = 5070, reps = 25))

base +
  geom_sf(data = eas[[2]][['Geometry']], aes(fill = factor(ID)), alpha = 0.2)

## ----ibd sampling, warning=F, message=F---------------------------------------
# create an arbitrary template for example - best to do this over your real range of all species collections
# so that it can be recycled across species. 
template <- terra::rast(terra::ext(bb), crs = terra::crs(spp), resolution = c(0.1, 0.1))
terra::values(template) <- 0

# the species range now gets 'burned' into the raster template. 
spp_concave <- purrr::map(spp_concave, \(x) {
  st_as_sf(x) |> 
    mutate(Range = 1, .before = geometry) |>
    st_transform(4326) |> 
    terra::rasterize(template, field = 'Range') 
})

## the actual sampling happens here within `map`
ibd_samples <- spp_concave |>
  purrr::map(~ IBDBasedSample(.x, n = 10, fixedClusters = FALSE, template = template, planar_proj = 5070))

base + ## visualize for a single taxon. 
  geom_sf(data = ibd_samples[[1]][['Geometry']], aes(fill = factor(ID)), alpha = 0.2)

## ----echo = F-----------------------------------------------------------------
rm(spp_concave)

## ----prioritize sample areas--------------------------------------------------
ibd_samples_priority <- ibd_samples |>
  purrr::map(~ st_transform(.x$Geometry, 5070)) |> 
  purrr::map(~ list(PrioritizeSample(.x, n_breaks = 3)))

base + 
  geom_sf(data = ibd_samples_priority[[2]][[1]][['Geometry']], aes(fill = factor(Level)), alpha = 0.2)

## ----remove bounding boxes, echo=FALSE----------------------------------------
rm(bb)

## ----write out data for storage, eval = F-------------------------------------
# ## create a directory to hold the outputs
# p2Collections <- file.path('~', 'Documents', 'WorkedExample_Output')
# dir.create(p2Collections, showWarnings = FALSE)
# 
# ## we will only save the template raster once since it is recycled across taxa.
# dir.create(file.path(p2Collections, 'IBD_raster_template'), showWarnings = FALSE)
# terra::writeRaster(template,
#   filename = file.path(p2Collections, 'IBD_raster_template', 'IBD_template.tif'), overwrite = FALSE)
# 
# ## save each species as a unique geopackage.
# for(i in seq_along(sppL)){
# 
#   fp = file.path(p2Collections, paste0(gsub(' ', '_', sppL[[i]]$species[1]), '.gpkg'))
# 
#   ### GBIF occurrence points
#   st_write(sppL[[i]], dsn = fp, layer =  'occurrence_points', quiet = TRUE)
# 
#   ### results of equal area sampling
#   st_write(eas[[i]]$Geometry,  dsn = fp, layer =  'equal_area_samples', quiet = TRUE, append = TRUE)
# 
#   ### ibd based sample (note can be reconstruced from the hulls of ibd samples priority)
#   st_write(ibd_samples[[i]]$Geometry, dsn = fp, layer =  'ibd_samples', quiet = TRUE, append = TRUE)
# 
#   ### prioritized information for the IBD samples.
#   st_write(ibd_samples_priority[[i]][[1]]$Geometry, dsn = fp, layer =  'ibd_sampling_priority', quiet = TRUE, append = TRUE)
# 
#   message(format(object.size(fp), standard = "IEC", units = "MiB", digits = 4))
# }

## ----clean up environment, echo = F, warning = F------------------------------
rm(spp, sppL, eas, ibd_samples, ibd_samples_priority, base, template, p2Collections, fp, i)

