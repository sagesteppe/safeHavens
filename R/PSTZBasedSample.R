#' Determine the sample size for n provisional seed transfer zones in an area
#' 
#' @description Intersect a vector data file of provisional seed transfer zones (PSTZs) (`pstz`) with the range of a focal taxon (`x`) and selects `n` seed zones (sz) to sample.
#' If fewer than n seed zones exist in the range of the species, then extra samples are added to the seed zones to meet the required `n` using a user specified method. 
#' 
#' 
#' This function has built in support for the official PSTZS FROM XXX created by the USFS, will work on other data sets with minimal user transformation to the input data and column specifications.  
#' If more seed zones exist across the range of the species than n (uncommon),  select a method to determine how samples are selected. 
#' Methods include: 
#'    - 1) `Largest` sampling from the largest seed zones by area
#'    - 2) `Smallest` smallest seed zones by area 
#'    - 3) `Most` the seed zones with the most disconnected habitat (measured by the number of polygons). 
#'    - 4) `Assist-warm` the seed zones with warmest climate 
#'    - 5) `Assist-drier`  the seed zones with the largest difference in warm and cold metrics
#' 
#' Each of these methods then returns the largest polygon within the criterion. 
#' 
#' A note on 'polygons'. 
#' Simple features are able to store polygon data in two main formats, a 'MULTIPOLYGON', where all individual polygons composing a class are stored collectively, or as 'POLYGONS' where each individual polygon is a unique entry within the class. 
#' 'Polygons' are generally used when two areas of the same class are discontinuous, and an analyst wants to easily analyze them separately.
#' 'MULTIPOLYGONS' are generally created by an analyst interested in understanding properties of the entire class, e.g. a single provisional seed zone. 
#' The function will test and cast 'MULTIPOLYGONS' to polygons, but this will incur overhead while iterating across data sets, it is best for you to cast the data to POLYGONS (as shown below).
#' 
#' @param x sf/tibble/data frame. A range of species as a simple feature (sf) object. 
#' @param n Numeric. desired total number of samples across this range
#' @param pstz sf/tibble/data frame. Provisional seed transfer zone vector data file (~shapefile) from PUBLICATION, loaded as an sf object. 
#' @param decrease_method Character. One of 'Largest', 'Smallest', 'Most', 'Assist-warm', 'Assist-drier'. See description for details. Default = 'Largest'. 
#' @param increase_method Character. One of 'Largest', 'Smallest', 'Most', 'Assist-warm', 'Assist-drier'. See description for details. Default = 'Largest'. 
#' @param warmest_col Character. Name of the variable with the information regarding WARMEST ... 
#' @param precip_range Character. Name of the variable with the information regarding COLEST ... 
#' @param group_col Character. Name of the grouping variable which identifies the unique pstzs. 
#' @examples \dontrun{
#' library(tidyverse)
#'
#' # define a square species-range
#' 
#' set.seed(23)
#' 
#' sr_mat <- rbind(
#'   c(0,0), c(10,0), c(10,10), c(0,10), c(0,0)
#' )
#' sr_poly <- sf::st_polygon(list(sr_mat))
#' x <- sf::st_sf(id = 1, geometry = sf::st_sfc(sr_poly))
#' 
#' rm(sr_mat, sr_poly)
#' 
#' pstz = data.frame(
#'   x = runif(10, min = -2, max = 12),
#'   y = runif(10, min = -2, max = 12)
#' ) |>
#'   sf::st_as_sf(coords = c('x', 'y')) |>
#'   sf::st_union() |>
#'   sf::st_voronoi() |>
#'   sf::st_collection_extract('POLYGON') |>
#'   sf::st_as_sf() |>
#'   dplyr::mutate(pstz_key = sample(LETTERS[1:7], size = 10, replace = T)) |>
#'   dplyr::rename('geometry' = x) |>
#'   sf::st_crop(x)
#' 
#' bp <- ggplot2::ggplot(x) + 
#'   geom_sf(fill = NA, lwd = 2) + 
#'   geom_sf(data = pstz, aes(fill = pstz_key)) 
#' 
#' bp + 
#'   geom_sf_label(data = pstz,aes(label = pstz_key))
#' 
#' 
#' # Example call #1: request number of samples == number 
#' 
#' ###################################################################### --- SUCCESS 
#' # example #1: request same numer of samples as zones - all zones returned. 
#' res1 <- sample_pstz(
#'   x = x, 
#'   n = length(unique(pstz[['pstz_key']])), 
#'   pstz = pstz, pstz_key = "pstz_key",
#'   increase_method = "Most"
#' )
#' 
#' bp +
#'   geom_sf(data = res1, alpha = 0.9) + 
#'   geom_sf_label(data = pstz,aes(label = pstz_key)) 
#' 
#' ## note that we get the largest polygon from EACH group to sample from. 
#' 
#' ##################################################################### --- SUCCESS
#' # Example #2: request fewer samples than zones -> subset by method - choosing largest by area 
#' 
#' res2 <- sample_pstz(x = x, n = 3, pstz = pstz, pstz_key = "pstz_key", increase_method = "Largest")
#' 
#' res2 |>
#'   group_by(pstz_key) |>
#'   mutate(total_area = sum(poly_area)) |>
#'   sf::st_drop_geometry() |>
#'   arrange(-total_area)|>
#'   knitr::kable()
#' 
#' bp + # picks, the two largest ## BUT NEED OT INVERT THE MASK 
#'   geom_sf(data = res2, alpha = 0.9) + 
#'   geom_sf_label(data = pstz,aes(label = pstz_key))
#' 
#' 
#' ####################################################################### -- SUCCESS 
#' # Example #3: request fewer samples than zones -> subset by method - choosing smallest by area 
#' res3 <- sample_pstz(x = x, n = 3, pstz = pstz, pstz_key = "pstz_key", increase_method = "Smallest")
#' 
#' res3 |>
#'   group_by(pstz_key) |>
#'   mutate(total_area = sum(poly_area)) |>
#'   sf::st_drop_geometry() |>
#'   arrange(total_area) |>
#'   knitr::kable()
#' 
#' ## note that we return the largest polygon (poly_area) within the `pstz_key` group, ranked by (total_area)
#' 
#' bp + # picks, the n smallest
#'   geom_sf(data = filter(res3, allocation == 0), alpha = 0.9) + 
#'   geom_sf_label(data = pstz, aes(label = pstz_key)) 
#' 
#' 
#' ####################################################################
#' # Example #4: request more samples than zones -> allocate extras to Largest pSTZs
#' 
#' ## note that is really a rounding rule - 'Largest' favors giving extra collections to the largest polygons
#' ## while 'smallest' favors giving them to smaller polygons. It is really mostly for edge cases, and
#' # the two will generally behave similarly on contrived examples. 
#' res4 <- sample_pstz(x = x, n = 12, pstz = pstz, pstz_key = "pstz_key", increase_method = "Largest")
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
#' 
#' ####################################################################
#' # Example #5: request more samples than zones -> allocate extras to pSTZs with most polygons
#' res5 <- sample_pstz(x = x, n = 14, pstz = pstz, pstz_key = "pstz_key", increase_method = "Most")
#'   
#' res5 |>
#'   group_by(pstz_key) |>
#'   summarize(Count = n(), Total_Allocation = sum(allocation)) |>
#'   sf::st_drop_geometry() |>
#'   arrange(-Count) |>
#'   knitr::kable()
#' 
#' }
#' @export
sample_pstz <- function(x, n = 20, pstz, pstz_key,
                        decrease_method = c("Largest","Smallest","Most","Assist-warm","Assist-drier"),
                        increase_method = c("Largest","Smallest","Most","Assist-warm","Assist-drier"),
                        warmest_col = NULL, precip_col = NULL) {
  decrease_method <- match.arg(decrease_method)
  increase_method <- match.arg(increase_method)
  
  if (missing(x) || missing(pstz) || missing(pstz_key)) {
    stop("Must supply x (range), pstz (zones), and pstz_key (group column name).")
  }

  # If climate-based methods are requested, check that the columns are provided.
  if (decrease_method %in% c("Assist-warm","Assist-drier") || 
      increase_method %in% c("Assist-warm","Assist-drier")) {
    if (decrease_method == "Assist-warm" || increase_method == "Assist-warm") {
      if (is.null(warmest_col) || !(warmest_col %in% names(pstz))) {
        stop("For Assist-warm, you must supply a valid warmest_col present in pstz.")
      }
    }
    if (decrease_method == "Assist-drier" || increase_method == "Assist-drier") {
      if (is.null(precip_col) || !(precip_col %in% names(pstz))) {
        stop("For Assist-drier, you must supply a valid precip_col present in pstz.")
      }
    }
  }
  
  sf::st_agr(x) <- "constant"
  sf::st_agr(pstz) <- "constant"
  
  # ensure all geometries are polygons (flatten MULTIPOLYGON or GEOMETRYCOLLECTIONs)
  pstz_poly <- sf::st_intersection(pstz, x)
  if (!all(sf::st_geometry_type(pstz_poly) == "POLYGON")) {
    pstz_poly <- sf::st_collection_extract(pstz_poly, "POLYGON", warn = FALSE)
  }
  pstz_poly <- sf::st_make_valid(pstz_poly)
  sf::st_agr(pstz_poly) <- "constant"

  # subset zones to those intersecting x
  pstz_sub <- sf::st_intersection(pstz_poly, x)
  if (!all(sf::st_geometry_type(pstz_sub) == "POLYGON")) {
    pstz_sub <- sf::st_collection_extract(pstz_sub, "POLYGON", warn = FALSE)
  }
  pstz_sub <- sf::st_make_valid(pstz_sub)
  
  if (nrow(pstz_sub) == 0) {
    warning("No PSTZ polygons intersect species range.")
    return(NULL)
  }
  
  # compute area for each polygon (units: e.g. m^2 by default; convert to hectares or keep as is)
  pstz_sub <- pstz_sub %>%
    dplyr::mutate(poly_area = units::set_units(sf::st_area(geometry), "m^2"))
  
  # compute zone-level summary: total area and polygon count
  zone_summary <- pstz_sub %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(!!sym(pstz_key)) %>%
    dplyr::summarise(
      Polygon_ct = n(),
      total_area_m2 = sum(poly_area), 
      .groups = 'drop'
    ) %>%
    dplyr::ungroup()
  
  # if climate-based methods, join the climate columns from original pstz to zone_summary
  # (assumes those columns are uniform per PSTZ key, or you define a rule)
  if ("Assist-warm" %in% c(decrease_method, increase_method)) {
    cs <- pstz_poly %>%
      sf::st_drop_geometry() %>%
      dplyr::select(!!sym(pstz_key), !!sym(warmest_col)) %>%
      dplyr::distinct()
    zone_summary <- zone_summary %>%
      dplyr::left_join(cs, by = pstz_key)
  }
  if ("Assist-drier" %in% c(decrease_method, increase_method)) {
    ps <- pstz_poly %>%
      sf::st_drop_geometry() %>%
      dplyr::select(!!sym(pstz_key), !!sym(precip_col)) %>%
      dplyr::distinct()
    zone_summary <- zone_summary %>%
      dplyr::left_join(ps, by = pstz_key)
  }

  # which zones exist in the intersection
  zones <- zone_summary[[pstz_key]]
  n_zones <- length(zones)
  
  if (n_zones == 0) {
    warning("No seed zones overlapping range after intersection.")
    return(NULL)
  }
  
  if (n == n_zones) {

    # one sample per zone  ; pick the largest polygon in each zone
    result <- pstz_sub %>%
      dplyr::group_by(!!sym(pstz_key)) %>%
      dplyr::slice_max(order_by = poly_area, n = 1) %>%
      dplyr::ungroup() |> 
      dplyr::mutate(allocation = 1, .before = 'geometry')

    return(result)
  }
  
  if (n < n_zones) {
    
    ranked <- rank_zones(zone_summary, decrease_method)
    selected <- ranked[1:n, ][[pstz_key]]
    
    result <- pstz_sub %>%
      dplyr::filter(.data[[pstz_key]] %in% selected) %>%
      dplyr::group_by(!!sym(pstz_key)) %>%
      dplyr::slice_max(order_by = poly_area, n = 1) %>%
      dplyr::mutate(allocation = 1, .before = 'geometry') |> 
      dplyr::ungroup()
    
    return(result)
  }  else  { # n > n_zones
  
    # case: n > n_zones → assign multiple samples per zone
    # strategy: assign one per zone, then assign the remainder to zones in descending polygon-size
    result <- allocate_increase(
      zone_summary = zone_summary,
      pstz_sub = pstz_sub,
      pstz_key = pstz_key,
      requested = n,
      method = increase_method,
      warmest_col = warmest_col,
      precip_col = precip_col
    )
    return(result)
    }
}

# Helper to rank zones based on method
#' @keywords internal
#' @noRd
rank_zones <- function(df, method) {
  # df must contain total_area_m2, Polygon_ct, plus optionally climate columns
  if (method == "Largest") {
    dplyr::arrange(df, desc(total_area_m2))
  } else if (method == "Smallest") {
    dplyr::arrange(df, total_area_m2)
  } else if (method == "Most") {
    dplyr::arrange(df, desc(Polygon_ct))
  } else if (method == "Assist-warm") {
    dplyr::arrange(df, desc(.data[[warmest_col]]))
  } else if (method == "Assist-drier") {
    df %>% dplyr::arrange(.data[[precip_col]])  # assume smaller precip → drier
  } else {
    stop("Unknown method: ", method)
  }
}

#' @keywords internal
#' @noRd
allocate_increase <- function(zone_summary, pstz_sub, pstz_key,
                              requested, method,
                              warmest_col = NULL, precip_col = NULL) {

  # ensure matching types
  if (!identical(class(pstz_sub[[pstz_key]]), class(zone_summary[[pstz_key]]))) {
    zone_summary[[pstz_key]] <- as.character(zone_summary[[pstz_key]])
    pstz_sub[[pstz_key]] <- as.character(pstz_sub[[pstz_key]])
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
      ranked <- dplyr::arrange(zone_summary, dplyr::desc(Polygon_ct))
      allocation_counts <- allocation_counts + rep(0, n_zones)  # base 1 + distribute remainder by rank
      allocation_counts[1:remainder] <- allocation_counts[1:remainder] + 1

    } else if (method == "Assist-warm") {
      ranked <- dplyr::arrange(zone_summary, dplyr::desc(!!sym(warmest_col)))
      allocation_counts <- allocation_counts[order(zone_summary[[pstz_key]] %in% ranked[[pstz_key]])]
      allocation_counts[1:remainder] <- allocation_counts[1:remainder] + 1

    } else if (method == "Assist-drier") {
      ranked <- dplyr::arrange(zone_summary, !!sym(precip_col))
      allocation_counts[1:remainder] <- allocation_counts[1:remainder] + 1

    } else {
      stop("Unknown increase_method: ", method)
    }
  }

  allocation <- tibble(
    !!pstz_key := zone_summary[[pstz_key]],
    n_alloc = allocation_counts
  )

  # distribute points across polygons per zone
  out_list <- vector("list", n_zones)
  for (i in seq_len(n_zones)) {
    zone_i <- allocation[[pstz_key]][i]
    points_to_assign <- allocation$n_alloc[i]

    polys <- pstz_sub %>% dplyr::filter(.data[[pstz_key]] == zone_i)
    n_polys <- nrow(polys)
    if (n_polys == 0) next

    base_poly <- points_to_assign %/% n_polys
    remainder_poly <- points_to_assign %% n_polys
    poly_alloc <- base_poly + c(rep(1, remainder_poly), rep(0, n_polys - remainder_poly))

    polys <- polys %>% dplyr::mutate(allocation = poly_alloc)
    out_list[[i]] <- polys
  }

  result <- do.call(rbind, out_list)
  return(result)
}

#' split and extract the temperature values from Tmin and AHM columns
#' 
#' @description Programmed for the Bower provisional seed zone products, a helper function for
#' separating and recovering the values from the columns in the data. Might be helpful on other data sets.  
#' 
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
    
  } else { df <- data.frame(t(ob))} # else transform straight to df

  df = setNames( # convert from characters to numbers for math
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
    median = apply(df, MARGIN = 1, FUN = function(x){median(x, na.rm = TRUE)}), 
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