setwd('~/Documents/assoRted/StrategizingGermplasmCollections/scripts')

polygon <- spData::us_states |>
  dplyr::select(NAME) |>
  dplyr::filter(NAME %in% c('California', 'Oregon')) |>
  sf::st_transform(4326)

ecoregions <- sf::st_read('../data/spatial/us_eco_l4/us_eco_l4_no_st.shp', quiet = TRUE) |>
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