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
#' @examples 
#' 
#' 
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
  
  # Helper to rank zones based on method
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
  
  if (n == n_zones) {
    # one sample per zone  ; pick the largest polygon in each zone
    result <- pstz_sub %>%
      dplyr::group_by(!!sym(pstz_key)) %>%
      dplyr::slice_max(order_by = poly_area, n = 1) %>%
      dplyr::ungroup()
    return(result)
  }
  
  if (n < n_zones) {
    # need to select a subset of zones (n of them), then one polygon per selected zone
    ranked <- rank_zones(zone_summary, decrease_method)
    selected <- ranked[1:n, ][[pstz_key]]
    
    result <- pstz_sub %>%
      dplyr::filter(.data[[pstz_key]] %in% selected) %>%
      dplyr::group_by(!!sym(pstz_key)) %>%
      dplyr::slice_max(order_by = poly_area, n = 1) %>%
      dplyr::ungroup()
    
    return(result)
  }
  
  # case: n > n_zones → assign multiple samples per zone
  # strategy: assign one per zone, then assign the remainder to zones in descending polygon-size
  base <- rep(1, n_zones)
  remainder <- n - n_zones
  
  ranked <- rank_zones(zone_summary, increase_method)
  # distribute extra: give +1 to top zones until remainder used up
  extra <- rep(0, n_zones)
  extra[1:remainder] <- 1
  allocation <- tibble(
    !!pstz_key := ranked[[pstz_key]],
    n_alloc = base + extra
  )
  
  # now for each zone, sample n_alloc polygons, selecting largest first (or fewer if zone has fewer)
  out_list <- list()
  for (i in seq_len(nrow(allocation))) {
    zone_i <- allocation[[pstz_key]][i]
    how_many <- allocation$n_alloc[i]
    polys <- pstz_sub %>% filter(.data[[pstz_key]] == zone_i)
    # if fewer polygons than allocate, just take them all (or sample duplicates? probably not)
    take_n <- min(nrow(polys), how_many)
    chosen <- polys %>%
      dplyr::slice_max(order_by = poly_area, n = take_n)
    out_list[[i]] <- chosen
  }
  
  result <- do.call(rbind, out_list)
  return(result)
}



library(tidyverse)

# define a square species-range

set.seed(23)

sr_mat <- rbind(
  c(0,0), c(10,0), c(10,10), c(0,10), c(0,0)
)
sr_poly <- sf::st_polygon(list(sr_mat))
x <- sf::st_sf(id = 1, geometry = sf::st_sfc(sr_poly))

rm(sr_mat, sr_poly)

pstz = data.frame(
  x = runif(10, min = -2, max = 12),
  y = runif(10, min = -2, max = 12)
) |>
  sf::st_as_sf(coords = c('x', 'y')) |>
  sf::st_union() |>
  sf::st_voronoi() |>
  sf::st_collection_extract('POLYGON') |>
  sf::st_as_sf() |>
  dplyr::mutate(pstz_key = sample(LETTERS[1:7], size = 10, replace = T)) |>
  dplyr::rename('geometry' = x) |>
  sf::st_crop(x)

bp <- ggplot2::ggplot(x) + 
  geom_sf(fill = NA, lwd = 2) + 
  geom_sf(data = pstz, aes(fill = pstz_key))

bp

# Example call #1: request number of samples == number of zones overlapping -> one per zone
res1 <- sample_pstz(x = x, n = 10, pstz = pstz, pstz_key = "pstz_key",
                    decrease_method = "Largest", increase_method = "Largest")
print(res1) 

bp + #  should pick the n largest zones overlapping species range
  geom_sf(data = res1, alpha = 0.9)

# Example call #2: request fewer samples than zones -> subset zones by method
res2 <- sample_pstz(x = x, n = 5, pstz = pstz, pstz_key = "pstz_key",
                    decrease_method = "Smallest", increase_method = "Largest")
print(res2)  # picks the smallest overlapping zone

bp + # picks, the two smallest zones. 
  geom_sf(data = res2, alpha = 0.9)

# Example call #3: request more samples than zones -> allocate extras
res3 <- sample_pstz(x = x, n = 15, pstz = pstz, pstz_key = "pstz_key",
                    decrease_method = "Largest", increase_method = "Largest")
print(res3)  










clim_df <- split_cols(sf::st_drop_geometry(pstz), 'Tmin_class', sep = '-') %>%
  dplyr::rename(Tmin_lower = lower, Tmin_upper = upper,
                Tmin_median = median, Tmin_range = range)

# (Optionally do the same for moisture / AHM / precip column.)

# Then bind to pstz (or join by pstz_key):
pstz2 <- pstz %>%
  dplyr::bind_cols(clim_df)

# Now pass pstz2 into sample_pstz(), with warmest_col = "Tmin_median" (or whatever you choose)
out <- sample_pstz(x = my_species_range, n = 20, pstz = pstz2,
                   pstz_key = "my_zone_id",
                   decrease_method = "Assist-warm",
                   warmest_col = "Tmin_median")




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

temps <- unlist_cols(df, 'Tmin_class')
moisure <- unlist_cols(df, 'AHM_class')



## sort by temperature, or moisture, 
order ( temps[['upper']] ) ## ascending , i.e. coolest to warmest
order ( temps[['upper']], decreasing = TRUE) ## warmest to coolest. 
 


## method for splitting points, when requested > zones. each stz gets one point + remainder
# this will result in the most even split between the number of points requested and the number of seed zones. 
requested = 20
n_zones <- 8
collections <- rep(0, n_zones)

base_amount <- requested %/% n_zones  # 1 divide, each elements get at least these. 
remainder <- requested %% n_zones      # these many accensions remain to be allocated. 
collections <- base_amount + # all get the base
  c(rep(1, remainder), rep(0, n_zones - remainder)) # the first few records get + 1, the others get no additional points


## now the initial input data is sorted. 




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
  
  proportions <- values / sum(values)
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








## increase method - less points requested than regions exist... 
requested = 4 
v = c(0.4, 0.3, 0.2, 0.1, 0.1)
length(v)

## option 1, select from the largest to smallest, as possible 
v [ order(v, decreasing = T)[seq_len(requested)] ]

## option 3, select from the smallest to largest, as possible
  v [ order(v) [seq_len(requested)] ]