setwd('~/Documents/assoRted/StrategizingGermplasmCollections/scripts')

library(tidyverse)
library(sf)

################################################################################
# First we do this in USA, so people have a standard Omernik L4 ecoregion shapefile
# to play with. This is the default for the function to utilize. 

polygon <- spData::us_states |>
  dplyr::select(NAME) |>
  dplyr::filter(NAME %in% c('California', 'Oregon')) |>
  sf::st_transform(4326)

paste0(colnames(ecoregions), collapse = ', ')

ecoregions <- sf::st_read('../data/spatial/us_eco_l4/us_eco_l4_no_st.shp', quiet = TRUE)# |>
  sf::st_transform(4326) |>
  sf::st_make_valid() |>
  sf::st_intersection(polygon, ) |>
  sf::st_cast('MULTIPOLYGON') |>
  rmapshaper::ms_simplify(keep = 0.01) |>
  sf::st_make_valid()

# OG file 95 MiB. -> 0.6 MiB
format(object.size(ecoregions), units = 'MiB')

setwd('~/Documents/assoRted/safeHavens/data/gpkg')
sf::st_write(ecoregions, 'WesternEcoregions.gpkg')

rm(ecoregions, polygon)


###############################################################################
# Now we do this for the Bradypus example for the vignette. 

neo_eco <- sf::st_read('../data/spatial/NeoTropics-Ecoregions/NeotropicMap_Geo.shp', quiet = TRUE) |>
  dplyr::filter(Provincias != 'N/A') |>
  sf::st_transform(4326) |>
  sf::st_make_valid() |>
  sf::st_cast('MULTIPOLYGON') |>
  rmapshaper::ms_simplify(keep = 0.01) |>
  sf::st_make_valid() 

format(object.size(neo_eco), units = 'MiB')
head(neo_eco)

setwd('~/Documents/assoRted/safeHavens/data/gpkg')
sf::st_write(neo_eco, 'NeoTropicsEcoregions.gpkg')

ggplot() + 
  geom_sf(data = neo_eco, aes(fill = Provincias))


   