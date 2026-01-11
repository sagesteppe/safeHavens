#' Sample spatial zones within a species range
#' 
#' @description Intersect a vector data file of spatial zones (ecoregions, provisional 
#' seed transfer zones, or other spatial partitions) with the range of a focal taxon 
#' and select `n` zones to sample. If fewer than n zones exist, extra samples are 
#' allocated using the specified method. If more than n zones exist, zones are selected 
#' using the specified method.
#' 
#' @details
#' Simple features can store polygon data as 'MULTIPOLYGON' (all polygons of a class 
#' stored collectively) or 'POLYGON' (each individual polygon is a unique entry). 
#' This function will cast MULTIPOLYGONS to POLYGONS as needed, but pre-casting 
#' will improve performance.
#' 
#' Available methods for zone selection:
#' - **Largest**: Select zones by total area (descending)
#' - **Smallest**: Select zones by total area (ascending)  
#' - **Most**: Select zones with most polygons (highest fragmentation)
#' - **Assist-warm**: Select warmest zones (requires `warmest_col`)
#' - **Assist-drier**: Select driest zones (requires `precip_col`)
#'
#' @param x sf object. Species range as a simple feature.
#' @param zones sf object. Spatial zones vector data (ecoregions, PSTZs, etc.).
#' @param zone_key Character. Column name identifying unique zones (required).
#' @param n Numeric. Desired total number of samples. Default = 20.
#' @param decrease_method Character. Method when n < number of zones. One of: 
#'   "Largest", "Smallest", "Most", "Assist-warm", "Assist-drier". Default = "Largest".
#' @param increase_method Character. Method when n > number of zones. One of:
#'   "Largest", "Smallest", "Most", "Assist-warm", "Assist-drier". Default = "Largest".
#' @param warmest_col Character. Column name for warmest temperature metric (required 
#'   for "Assist-warm" method).
#' @param precip_col Character. Column name for precipitation metric (required for 
#'   "Assist-drier" method).
#'   
#' @returns sf object with selected zones/polygons and an `allocation` column 
#'   indicating number of samples per polygon.
#'   
#' @examples
#' \dontrun{
#' library(tidyverse)
#' 
#' sr_mat <- rbind(
#'   c(0,0), c(10,0), c(10,10), c(0,10), c(0,0)
#' )
#' sr_poly <- sf::st_polygon(list(sr_mat))
#' x <- sf::st_sf(id = 1, geometry = sf::st_sfc(sr_poly))
#' 
#' rm(sr_mat, sr_poly)
#' 
#' zone_polys = data.frame(
#'   ## randomly generate some points in XY space. 
#'   x = runif(10, min = -2, max = 12),
#'   y = runif(10, min = -2, max = 12)
#' ) |>
#'   # conver to spatial points
#'   sf::st_as_sf(coords = c('x', 'y')) |>
#'   ## allocate XY space to it's nearest point
#'   sf::st_union() |>
#'   sf::st_voronoi() |>
#'   ## extract the contiguous pieces of XY space around points
#'   sf::st_collection_extract('POLYGON') |>
#'   sf::st_as_sf() |>
#'   ## make up seed zones on the fly, assign multiple polygons to some zones. 
#'   dplyr::mutate(pstz_key = sample(LETTERS[1:7], size = 10, replace = T)) |>
#'   dplyr::rename('geometry' = x) |>
#'   sf::st_crop(x)
#' 
#' bp <- ggplot2::ggplot(x) + 
#'   ggplot2::geom_sf(fill = NA, lwd = 2) + 
#'   ggplot2::geom_sf(data = zone_polys, ggplot2::aes(fill = pstz_key)) 
#' 
#' bp + 
#'   ggplot2::geom_sf_label(data = zone_polys, ggplot2::aes(label = pstz_key))
#' 
#' ###################################################################### 
#' # example #1: request same numer of samples as zones - all zones returned. 
#'  res1 <- PolygonBasedSample(
#'    x = x, 
#'    n = length(unique(zone_polys[['pstz_key']])), 
#'    zones = zone_polys, 
#'    zone_key  = "pstz_key",
#'    increase_method = "Most"
#'  )
#' 
#' bp +
#'   geom_sf(data = res1, alpha = 0.9) + 
#'   geom_sf_label(data = pstz,aes(label = pstz_key)) 
#' 
#' ## note that we get the largest polygon from EACH group to sample from. 
#' 
#' ##################################################################### 
#' # Example #2: request fewer samples than zones -> subset by method - choosing largest by area 
#' 
#' res2 <- PolygonBasedSample(
#'    x = x, n = 3, zones = zone_polys, zone_key  = "pstz_key", increase_method = "Largest"
#' )
#' 
#' res2 |>
#'   group_by(pstz_key) |>
#'   mutate(total_area = sum(poly_area)) |>
#'   sf::st_drop_geometry() |>
#'   arrange(-total_area)|>
#'   knitr::kable()
#' 
#' bp + # picks, the three largest 
#'   geom_sf(data = res2, alpha = 0.9) + 
#'   geom_sf_label(data = zone_polys, aes(label = pstz_key))
#' 
#' ####################################################################### 
#' # Example #3: request fewer samples than zones -> subset by method - choosing smallest by area 
#' res3 <- PolygonBasedSample(
#'   x = x, n = 3, zones = zone_polys, zone_key = "pstz_key", increase_method = "Smallest")
#' 
#' res3 |>
#'   group_by(pstz_key) |>
#'   mutate(total_area = sum(poly_area)) |>
#'   sf::st_drop_geometry() |>
#'   arrange(total_area) |>
#'   knitr::kable()
#' 
#' ## returns the largest polygon (poly_area) within the `pstz_key` group, ranked by (total_area)
#' 
#' bp + # picks, the n smallest - too small to see sometimes
#'   geom_sf(data = filter(res3, allocation == 0), alpha = 0.9) + 
#'   geom_sf_label(data = zone_polys, aes(label = pstz_key)) 
#' 
#' ####################################################################
#' # Example #4: request more samples than zones -> allocate extras to Largest pSTZs
#'
#' ## note that is really a rounding rule - 'Largest' favors giving extra collections to the largest 
#' ## polygons while 'smallest' favors giving them smaller polygons. It is really mostly for edge cases
#' ## and the two will generally behave similarly on contrived examples. 
#' res4 <- PolygonBasedSample(
#'    x = x, n = 12, zones = zone_polys, zone_key = "pstz_key", increase_method = "Largest")
#'   
#' res4 |>
#'   group_by(pstz_key) |>
#'   summarize(total_area = sum(poly_area),  Total_Allocation = sum(allocation)) |>
#'   sf::st_drop_geometry() |>
#'   arrange(-total_area) |>
#'   knitr::kable()
#' 
#' bp + 
#'   theme(legend.position = 'none') + 
#'   geom_sf(data = res4, aes(fill = as.factor(allocation))) + 
#'   geom_sf_label(data = res4, aes(label = allocation)) 
#' 
#' ####################################################################
#' # Example #5: request more samples than zones -> allocate extras to pSTZs with most polygons
#' res5 <- PolygonBasedSample(
#'    x = x, n = 14, zones = zone_polys, zone_key = "pstz_key", increase_method = "Most")
#'   
#' res5 |>
#'   group_by(pstz_key) |>
#'   summarize(Count = n(), Total_Allocation = sum(allocation)) |>
#'   sf::st_drop_geometry() |>
#'   arrange(-Count) |>
#'   knitr::kable()
#' }
#' @export
PolygonBasedSample <- function(
    x, 
    zones, 
    zone_key,
    n = 20,
    decrease_method = c("Largest", "Smallest", "Most", "Assist-warm", "Assist-drier"),
    increase_method = c("Largest", "Smallest", "Most", "Assist-warm", "Assist-drier"),
    warmest_col = NULL,
    precip_col = NULL
) {
  
  # Input validation
  decrease_method <- match.arg(decrease_method)
  increase_method <- match.arg(increase_method)
  
  if (missing(x) || missing(zones) || missing(zone_key)) {
    stop("Must supply `x` (species range), `zones` (spatial data for Ecoregions or STZs), and `zone_key` (grouping column).")
  }
  
  if (!zone_key %in% names(zones)) {
    stop("zone_key `", zone_key, "` not found in `zones` data.")
  }
  
   # Ensure all geometries are polygons
   # quelch warnings from s2 engine. 
  if(inherits(x, 'sf')){sf::st_agr(x) <- "constant"}
  if(inherits(zones, 'sf')){sf::st_agr(zones) <- "constant"}

  # Validate climate columns for climate-based methods
  climate_methods <- c("Assist-warm", "Assist-drier")
  if (decrease_method %in% climate_methods || increase_method %in% climate_methods) {
    if (decrease_method == "Assist-warm" || increase_method == "Assist-warm") {
      if (is.null(warmest_col) || !warmest_col %in% names(zones)) {
        stop("Assist-warm requires valid `warmest_col` variable in `zones` data.")
      }
    }
    if (decrease_method == "Assist-drier" || increase_method == "Assist-drier") {
      if (is.null(precip_col) || !precip_col %in% names(zones)) {
        stop("Assist-drier requires valid `precip_col` variable in `zones` data.")
      }
    }
  }
  

  zones_poly <- zones
  if (!all(sf::st_geometry_type(zones_poly) == "POLYGON")) {
    zones_poly <- sf::st_collection_extract(zones_poly, "POLYGON", warn = FALSE)
  }
  zones_poly <- sf::st_make_valid(zones_poly)

  sf::st_agr(zones_poly) <- "constant"
  
  # Intersect zones with species range
  zones_sub <- sf::st_intersection(zones_poly, x)
  if (!all(sf::st_geometry_type(zones_sub) == "POLYGON")) {
    zones_sub <- sf::st_collection_extract(zones_sub, "POLYGON", warn = FALSE)
  }
  zones_sub <- sf::st_make_valid(zones_sub)
  
  if (nrow(zones_sub) == 0) {
    warning("No zone polygons intersect species range.")
    return(NULL)
  }
  
  # Compute polygon areas
  zones_sub <- dplyr::mutate(zones_sub, poly_area = as.numeric(sf::st_area(geometry))) # in m^2
  
  # Compute zone-level summaries
  zone_summary <- zones_sub |>
    sf::st_drop_geometry()  |>
    dplyr::group_by(.data[[zone_key]]) |>
    dplyr::summarise(
      polygon_ct = dplyr::n(),
      total_area_m2 = sum(poly_area),
      .groups = 'drop'
    )
    
  
  ## prep and add climate data as required. 
  need_warm  <- "Assist-warm"  %in% c(decrease_method, increase_method)
  need_drier <- "Assist-drier" %in% c(decrease_method, increase_method)
  if (need_warm || need_drier) {

    modifiers_reduced <- zones_sub |>
      sf::st_drop_geometry() |>
      dplyr::group_by(.data[[zone_key]]) |>
      dplyr::summarise(
        !!!c(
          if (need_warm) {
            rlang::set_names(
              list(
                rlang::expr(
                  max(.data[[!!warmest_col]], na.rm = TRUE)
                )
              ),
              warmest_col
            )
          },
          if (need_drier) {
            rlang::set_names(
              list(
                rlang::expr(
                  min(.data[[!!precip_col]], na.rm = TRUE)
                )
              ),
              precip_col
            )
          }
        ),
        .groups = "drop"
      )

  zone_summary <- dplyr::left_join(zone_summary, modifiers_reduced, by = zone_key)
  }
  
  # Get unique zones
  unique_zones <- zone_summary[[zone_key]]
  n_zones <- length(unique_zones)
  
  if (n_zones == 0) {
    warning("No zones overlapping range after intersection.")
    return(NULL)
  }
  
  # Case 1: n == n_zones (one sample per zone)
  if (n == n_zones) {
    result <- zones_sub |>
      dplyr::group_by(.data[[zone_key]]) |>
      dplyr::slice_max(order_by = poly_area, n = 1) |>
      dplyr::ungroup() |>
      dplyr::mutate(allocation = 1, .before = geometry)
    return(result)
  }
  
  # Case 2: n < n_zones (select subset of zones)
  if (n < n_zones) {
    ranked <- rank_zones(zone_summary, decrease_method, total_area_m2, polygon_ct, warmest_col, precip_col)
    selected <- ranked[1:n, ][[zone_key]]
    
    result <- zones_sub |>
      dplyr::filter(.data[[zone_key]] %in% selected) |>
      dplyr::group_by(.data[[zone_key]]) |>
      dplyr::slice_max(order_by = poly_area, n = 1) |>
      dplyr::ungroup() |>
      dplyr::mutate(allocation = 1, .before = geometry)
    
    return(result)
  }
  
  # Case 3: n > n_zones (allocate multiple samples per zone)
  result <- allocate_increase(
    zone_summary = zone_summary,
    zones_sub = zones_sub,
    zone_key = zone_key,
    requested = n,
    method = increase_method,
    warmest_col = warmest_col,
    precip_col = precip_col
  )
}


# Helper to rank zones based on method
#' @keywords internal
#' @noRd
rank_zones <- function(df, method, total_area_m2 = NULL, Polygon_ct = NULL, warmest_col = NULL, precip_col = NULL){
  # df must contain total_area_m2, polygon_ct, plus optionally climate columns
  if (method == "Largest") {
    dplyr::arrange(df, dplyr::desc(total_area_m2))
  } else if (method == "Smallest") {
    dplyr::arrange(df, total_area_m2)
  } else if (method == "Most") {
    dplyr::arrange(df, dplyr::desc(polygon_ct))
  } else if (method == "Assist-warm") {
    dplyr::arrange(df, dplyr::desc(.data[[warmest_col]]))
  } else if (method == "Assist-drier") {
    dplyr::arrange(df, .data[[precip_col]])  # assume lower precip = drier. not globally safe. 
  } else {
    stop("Unknown method: ", method)
  }
}

#' @keywords internal
#' @noRd
allocate_increase <- function(zone_summary, zones_sub, zone_key, requested, method, warmest_col = NULL, precip_col = NULL) {

  # ensure matching types
  if (!identical(class(zones_sub[[zone_key]]), class(zone_summary[[zone_key]]))) {
    zone_summary[[zone_key]] <- as.character(zone_summary[[zone_key]])
    zones_sub[[zone_key]] <- as.character(zones_sub[[zone_key]])
  }

  n_zones <- nrow(zone_summary)
  if (n_zones == 0) return(NULL)

  allocation_counts <- rep(1, n_zones)  # start with 1 per zone
  remainder <- requested - n_zones      # remaining points to distribute

  if (remainder > 0) {
    if (method %in% c("Largest", "Smallest")) {

      # compute proportional allocation based on total_area_m2
      allocation_counts <- allocation_counts +
        proportional_round(
          values = zone_summary$total_area_m2,
          target_sum = remainder,
          method = ifelse(method == "Largest", "larger_up", "larger_down")
        )
    } else if (method == "Most") {

      top_indices <- order(zone_summary$polygon_ct, decreasing = TRUE)[1:min(remainder, n_zones)]
      allocation_counts[top_indices] <- allocation_counts[top_indices] + 1

    } else if (method == "Assist-warm") {

        ord <- order(zone_summary[[warmest_col]], decreasing = TRUE)
        for (i in seq_len(remainder)) {
          idx <- ord[(i - 1) %% n_zones + 1]
          allocation_counts[idx] <- allocation_counts[idx] + 1
        }

    } else if (method == "Assist-drier") {

        ord <- order(zone_summary[[precip_col]], decreasing = FALSE)
        for (i in seq_len(remainder)) {
          idx <- ord[(i - 1) %% n_zones + 1]
          allocation_counts[idx] <- allocation_counts[idx] + 1
        }
      
    } else {
      stop("Unknown increase_method: ", method)
    }
  }

  allocation <- data.frame(
    n_alloc = allocation_counts,
    stringsAsFactors = FALSE
  )
  allocation[[zone_key]] <- zone_summary[[zone_key]]


  # distribute points across polygons per zone
  out_list <- vector("list", n_zones)
  for (i in seq_len(n_zones)) {
    zone_i <- allocation[[zone_key]][i]
    points_to_assign <- allocation$n_alloc[i]

    polys <- dplyr::filter(zones_sub, .data[[zone_key]] == zone_i)
    n_polys <- nrow(polys)
    if (n_polys == 0) next

    base_poly <- points_to_assign %/% n_polys
    remainder_poly <- points_to_assign %% n_polys
    poly_alloc <- base_poly + c(rep(1, remainder_poly), rep(0, n_polys - remainder_poly))

    polys <- dplyr::mutate(polys, allocation = poly_alloc)
    out_list[[i]] <- polys
  }

  result <- do.call(rbind, out_list)
}
  
#' split and extract the temperature values from Tmin and AHM columns
#' 
#' @description Programmed for the Bower provisional seed zone products, a helper function for
#' separating and recovering the values from the columns in the data.
#' @param dat data frame with the columns required to split.
#' @param y character. column name to split.
#' @param sep character. separator between the values to split on. Default is '-'.
#' @examples
#' df = data.frame(
#'   'Tmin_class' = c('10 - 15 Deg. F.', '15 - 20 Deg. F.', '> 55 Deg. F.' ),
#'   'AHM_class' = c('2 - 3', '6 - 12', '3 - 6')
#' )
#' split_cols(df, 'Tmin_class')
#' split_cols(df, 'AHM_class')
#' 
#' @export
## unlisting this and collapse is super gross in base. but swollen on packages. 
split_cols <- function(dat, y, sep = '-'){

  ob <- sapply(strsplit(dat[[y]], sep), \(x) gsub('\\D+','', x))

  if(inherits(ob, 'list')){ # if unequal length information use this. 

    df <- t(as.data.frame(do.call(cbind, lapply(ob, \(x) {
      c(x, rep(NA, max(sapply(ob, length)) - length(x)))
    }))))
    
  } else { df <- data.frame(t(ob))} # else transform straight to data frame

  df = stats::setNames( # convert from characters to numbers for math
    data.frame(apply(df, FUN = as.numeric, MARGIN = 2)), 
    nm = c('lower', 'upper')
  )

  ## if the upper value is missing, make it the same as the lower. 
  if(any(is.na(df['upper']))){ # impute the upper value - assume same as lower.  
    idx <- which(is.na(df['upper']))
    df[idx, 'upper'] <- df[idx, 'lower']
  }

  data.frame( ## return object
    df,
    median = apply(df, MARGIN = 1, FUN = function(x){stats::median(x, na.rm = TRUE)}), 
    range = df[,2] - df[,1]
  )
}

#' round values up or down, based on their relative size. 
#' @examples
#' v = c(0.4, 0.3, 0.2, 0.1)
#' v*15
#' proportional_round(v, 15, method = "larger_up")
#' proportional_round(v, 15, method = "larger_down")
#' 
#' @keywords internal
#' @noRd
proportional_round <- function(values, target_sum, method = c("larger_up", "larger_down")) {
  method <- match.arg(method)
  
  proportions <- as.numeric( values / sum(values))
  ideal <- proportions * target_sum
  
  if (method == "larger_up") {
    # Start with floor, round up largest remainders
    result <- floor(ideal)
    remainders <- ideal - result
    shortage <- target_sum - sum(result)
    
    if (shortage > 0) {
      indices_to_adjust <- order(remainders, decreasing = TRUE)[1:shortage]
      result[indices_to_adjust] <- result[indices_to_adjust] + 1
    }
    
  } else {  # "larger_down"
    # Start with ceiling, round down largest remainders
    result <- ceiling(ideal)
    remainders <- result - ideal
    excess <- sum(result) - target_sum
    
    if (excess > 0) {
      indices_to_adjust <- order(remainders, decreasing = TRUE)[1:excess]
      result[indices_to_adjust] <- result[indices_to_adjust] - 1
    }
  }
  return(result)
}