#' Determine the sample size for n ecoregions in an area
#' 
#' @description Intersect a vector data file of ecoregions (`ecoregions`) with the range of a focal taxon (`x`) and selects `n` ecoregions to sample.
#' If fewer than n ecoregions exist in the range of the species, then extra samples are added to the largest ecoregions by their
#' area proportions to meet n. 
#' 
#' This function has adequate support for the official Omernik L4 shapefiles created by the EPA, and has worked with minimal experimentation to other vector data sources. 
#' If more ecoregions exist across the range of the species than n (very common),  select a method to determine how sample sizes (in this instance 1 sample per ecoregion) are selected. 
#' Methods include sampling from the n largest, orsmallest ecoregions by number, or the ecoregion with the most disconnected habitat (measured soley by the number of polygons). 
#' Each of these three methods then returns the largest polygon within the criterion. 
#'  
#'  A note on 'polygons'. 
#'  Simple features are able to store polygon data in two main formats, a 'MULTIPOLYGON', where all individual polygons composing a class are stored collectively, or as 'POLYGONS' where each individual polygon is a unique entry within the class. 
#'  'Polygons' are generally used when two areas of the same class are discontinuous, 
#'  and an analyst wants to easily analyze them separately.
#'  'MULTIPOLYGONS' are generally created by an analyst interested in understanding properties of the entire class. 
#'  The EPA Omernik spatial data set comes with both 'POLYGONS' and 'MULTIPOLYGONS', I have used it somewhat extensively and believe that the creators struck a happy balance between creating too many small polygons, e.g.  for areas like coastal reef island (a MULTIPOLYGON use case), and big polygons. 
#'  I do not modify them here, and on rare occasion (essentially islands), what I refer to as a 'polygon' may technically be a multipolygon. 
#'  
#' @param x a range of species as a simple feature (sf) object. 
#' @param ecoregions An ecoregion vector data file (~shapefile) in sf. 
#' @param OmernikEPA Boolean. TRUE indicates that the data are from the US EPA and minimally modified, if FALSE several more input are required to ensure the function maps over appropriately. 
#' If left blank the default method is to scan for an exact match of the standard Omernik ecoregion field (column) names, and if matched dispatched to the Omernik module, else fail unless the associated columns are specified (SEE BELOW). 
#' @param ecoregion_col Character. 
#' Name of the column contain the finest resolution data which are to be used for the analysis.
#'  For an Omernik L4 file it defaults to the relevant columns automatically, if a different type of file is submitted, and these not specified the function fails. 
#' @param n Numeric. desired total number of samples across this range
#' @param increase_method Character. Method to implement if the number of L4 ecoregions is
#' less than n. 
#' @param decrease_method Character. Method to implement if the number of L4 ecoregions is
#' greater than n. 
#' 'Largest' (the default) will select the n largest ecoregions by total area, and then select the largest single polygon within each of these classes. 
#' 'Smallest' will select the n smallest ecoregions by total area, and then select the largest single polygon within these classes. 
#' 'Most' will select the n ecoregions with the most polygons, and select the largest polygon from each. 
#' @examples 
#' 
#' # First example is using a subset (and with simplified geometry) Omernik L4
#' # ecoregion shapefile from the EPA. Which as of the time of writing were 
#' # available, at no cost, at the following URL
#' # https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states
#'
#' polygon <- spData::us_states |>
#' dplyr::select(NAME) |>
#'    dplyr::filter(NAME == 'California') |>
#'    sf::st_transform(4326)
#'    
#' Weco <- sf::st_read(
#'   file.path(system.file(package="safeHavens"), 'extdata', 'WesternEcoregions.gpkg'), 
#'   quiet = TRUE)
#' head(Weco)
#' 
#' out <- EcoregionBasedSample(polygon, Weco)
#' sum(out$n)
#' 
#' ggplot2::ggplot() + 
#'    ggplot2::geom_sf(data = out, ggplot2::aes(fill = factor(n)))
#'   
#' # This second example is from a recent publication by Morreno et al. 2022 and 
#' # presents biogeographic regions of the Neotropics and is available from a 
#' # google drive linked in the publication describing there creation located at
#' # https://www.scielo.br/j/aabc/a/hPft4CK6RV8QBr8nP7bxhRQ/?lang=en#
#' 
#' # Essentially we showcase how a user can maintain this functions utility
#' # while catering to data in a format differing from the Omernik L4 distribution. 
#' 
#' neo_eco <- sf::st_read(
#'    file.path(system.file(package="safeHavens"), 'extdata', 'NeoTropicsEcoregions.gpkg'),
#'     quiet = TRUE)
#'sp_range <- sf::st_polygon( # complete
#'  list(
#'    rbind(
#'      c(-80,-5), c(-80,10), c(-60,10), c(-55,5),
#'      c(-60,-5), c(-80,-5) 
#'    )
#'  )
#' ) |>
#'  sf::st_sfc() |> 
#'  sf::st_as_sf() |>
#'  sf::st_set_crs(4326) |>
#'  dplyr::rename(geometry = x) |>
#'  dplyr::mutate(Species = 'Da species')
#'
#' out <- EcoregionBasedSample(sp_range, neo_eco, ecoregion_col = 'Provincias')
#' sum(out$n)
#'
#'ggplot2::ggplot() + 
#'  ggplot2::geom_sf(data = neo_eco) + 
#'  ggplot2::geom_sf(data = out, ggplot2::aes(fill = factor(n))) + 
#'  ggplot2::geom_sf(data = sp_range, fill = NA, color = 'Red') 
#' 
#' # Note that both of the files of the above ecoregions have had their geometry
#' # simplified, i.e. made less complex - you should notice they look slightly angular
#' # like an old cartoon such as Rugrats or so. We do this to reduce the file size
#' # to make it easier to install the package, and reduce the run time of the functions
#' # for these simple examples. 
#' 
#' @returns An sf object, the same length as the input data set, with only the finest resolution eco level, and geometry fields retained, and a new column 'n' indicating how many accession should be gather from the ecoregion. 
#' @export 
EcoregionBasedSample <- function(x, ecoregions, OmernikEPA, n, ecoregion_col, increase_method, decrease_method){
  
  if(missing(OmernikEPA) & missing(ecoregion_col)){
    if(any(colnames(ecoregions) %in% 'L4_KEY')){OmernikEPA <- TRUE} else {
      OmernikEPA <- FALSE
      stop('A difference between the user submitted ecoregions data and standard Omernik L4 ecoregions has been detected. Please ensure you provide the appropriate colnames to fns `ecoregion_col`, set `OmernikEPA` to FALSE and try again.')}
  }
  
  if(missing(ecoregion_col)){L4_KEY <- 'L4_KEY'} else {L4_KEY <- ecoregion_col}
#  L4_KEY_quo <- rlang::enquo(L4_KEY)
  if(missing(n)){n<-20} 
  if(missing(increase_method)){increase_method<-'Area'}
  if(missing(decrease_method)){decrease_method<-'Largest'}
  sf::st_agr(ecoregions) = 'constant'
  
  # reduce their input ecoregions data set to the areas where x is. 
  
  ecoregions <- ecoregions  |> 
    sf::st_make_valid() |> 
    dplyr::mutate(ID = 1:dplyr::n(), .before = 1) 
  
  ecoregions_sub <- sf::st_intersection(ecoregions, x) |>
    sf::st_make_valid()
  ecoregions_sub <- ecoregions_sub [ !sf::st_is_empty(ecoregions_sub), ]
  ecoregions_sub <- ecoregions_sub[sf::st_geometry_type(ecoregions_sub) %in% c('POLYGON', 'MULTIPOLYGON'),]
  
  area <- dplyr::mutate(
    ecoregions_sub, Area = units::set_units(sf::st_area(ecoregions_sub), ha), .before = dplyr::last_col())
  
  area_summaries <- area |> 
    sf::st_drop_geometry() |> 
    dplyr::group_by(!!rlang::sym(L4_KEY)) |> 
    dplyr::summarise(
      Eco_lvl = 4,
      Polygon_ct = dplyr::n(), 
      Total_area = sum(Area)
    ) |>
    dplyr::rename(Name = !!rlang::sym(L4_KEY)) |>
    dplyr::ungroup()
  
  eco_lvls_ct <- data.frame(
    Eco_lvl = 'L4', 
    ct =  length(unique(area[,L4_KEY]))
  )
  
  cols <- c('ID', L4_KEY, 'n', 'geometry', 'x')
  
  # in this method we only assign counts to each area. 
  if(eco_lvls_ct[eco_lvls_ct$Eco_lvl=='L4','ct'] == n){
    
    if(sum(area_summaries[area_summaries, 'Polygon_ct']) == n){
      # if n polygons == n  woohoo! each ecoregion is allocated a sample size of one point 
      
      out <- dplyr::mutate(polygons, n = 1) |>
        dplyr::select(dplyr::any_of(cols))
      
    } else {
      # if n polygons > n, select a polygon in each ecoregion to represent it.  Either by AREA or CENTRALITY 
      if(increase_method == 'Area'){
        
        out <- dplyr::group_by(polygons, !!L4_KEY_quo) |>
          dplyr::arrange(Area, .by_group = TRUE) |>
          dplyr::slice_max(n = 1) |>
          dplyr::mutate(n = 1) |>
          dplyr::select(dplyr::any_of(cols))
        
      } else {
        
        unions <- dplyr::group_by(polygons, !!rlang::sym(L4_KEY)) |>
          dplyr::summarise(geometry = sf::st_union(geometry))
        pts <- sf::st_point_on_surface(unions) 
        out <- polygons[lengths(sf::st_intersects(pts, unions))>0, ] |>
          dplyr::mutate(n = 1) |>
          dplyr::select(dplyr::any_of(cols))
        
      }
    }
    
  } else if(eco_lvls_ct[eco_lvls_ct$Eco_lvl=='L4','ct'] < n) {
    # each L4 ecoregion is represented by a polygon, fewer than n L4's also means 
    # fewer than n polygons. The only option in this scenario is to add multiple
    # points per polygon based on AREA. 
    
    pct_area <- as.numeric(area$Area / sum(area$Area)) * 100
    sample <- numeric(length = length(pct_area))
    sample[pct_area<5] <- 1 # these small areas by default get a point. 
    n_remain <- n - sum(sample)
    
    if(n_remain >= n/2){
      pct_record <- sum(pct_area[pct_area>5]) / n_remain
      sample[pct_area>5] <- round(pct_area[pct_area>5] / pct_record)
      out <- dplyr::mutate(area, n = sample) |>
        dplyr::select(dplyr::any_of(cols))
    } else {
      # some areas bay be composed entirely of very small coverage areas. 
      # the largest 20 will get the records in that case. 
      area$n <- NA
      area$n[(order(pct_area, decreasing = TRUE)[1:20])] <- 1
      out <- dplyr::select(area, dplyr::any_of(cols))
    }
    
  } else { # MANY L4's across the species range, 
    
    # offer three options for how to select the target L4's
    # 1 & 2, by area, either descending - the largest ones, or ascending - the smallest ones
    # by the number of unique polygons per L4, which can be roughly considered as
    # representing the amount of discontinuity of each L4. 
    
    if(decrease_method=='Largest'){
      out <- area[area$L4_KEY %in% area_summaries[order(area_summaries$Total_area, 
                                                        decreasing = TRUE)[1:n],]$Name, ]
      out <- dplyr::group_by(out, L4_KEY) |> 
        dplyr::arrange(dplyr::desc(Area), .by_group = TRUE) |>
        dplyr::slice(1) 
    } else if(decrease_method=='Smallest'){
      out <- area[area$L4_KEY %in% area_summaries[order(area_summaries$Total_area,
                                                        decreasing = FALSE)[1:n],]$Name, ]
      out <- dplyr::group_by(out, L4_KEY) |> 
        dplyr::arrange(dplyr::desc(Area), .by_group = TRUE) |>
        dplyr::slice(1) 
    } else {
      out <- area[area$L4_KEY %in% area_summaries[order(area_summaries$Polygon_ct,
                                                        decreasing = TRUE)[1:n],]$Name, ]
      out <- dplyr::group_by(out, L4_KEY) |> 
        dplyr::arrange( dplyr::desc(Area), .by_group = TRUE) |>
        dplyr::slice(1) 
    }
    
    out <- dplyr::mutate(out, n = 1) |>
      dplyr::select(dplyr::any_of(cols)) 
  }
  
  # above we only added the desired sample sizes to each targeted ecoregion, 
  # now we will add back the ecoregions which are not targeted for sampling to reach
  # the n = 20 goal. 
  cols <- c('ID', L4_KEY, 'n', 'geometry')
  
  out_assigned <- out[!is.na(out$n),]
  ids <- unique(dplyr::pull(out_assigned, ID))
  
  out <- ecoregions |>
    dplyr::mutate(n = dplyr::if_else(ID %in% ids, 1, 0)) |>
    dplyr::arrange(ID) |>
    dplyr::select(dplyr::any_of(cols)) |>
    dplyr::select(-dplyr::any_of('ID')) |>
    sf::st_make_valid() |>
    dplyr::rename(dplyr::any_of( c(geometry = 'geom')))
  
  out <- sf::st_intersection(out, x)
  
  return(out)
}
