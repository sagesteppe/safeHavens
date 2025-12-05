# Example usage
set.seed(20)
library(ggplot2)

n_sites <- 50 # number of known populations
df <- data.frame(
  site_id = seq_len(n_sites),
  lat = runif(n_sites, 25, 30), # play with these to see elongated results. 
  lon = runif(n_sites, -125, -120),
  required = FALSE,
  coord_uncertainty = 0
)

## we will simulate coordinate uncertainty on a number of sites.  
uncertain_sites <- sample(setdiff(seq_len(n_sites), which(df$required)), size = min(6, n_sites-3))
df$coord_uncertainty[uncertain_sites] <- runif(length(uncertain_sites), 1e4, 1e5) # meters

## assign a required point. 
dists2c <- greatCircleDistance(
  median(df$lat), 
  median(df$lon), 
  df$lat, 
  df$lon
)
df[order(dists2c)[1],'required'] <- TRUE

dist_mat <- sapply(1:nrow(df), function(i) {
  greatCircleDistance(
    df$lat[i], df$lon[i],
    df$lat, df$lon
  )
})

dat = list(
  distances = as.matrix(dist_mat),
  sites = df
  )

resin <- maximizeDispersion(dat, n_sites = 5,  n_bootstrap=99, dropout_prob = 0.1)
resin$stability_score

ggplot(data = resin$input_data, 
  aes(
    x = lon, 
    y = lat, 
    shape = required, 
    size = cooccur_strength,
    color = selected
    )
  ) +
  geom_point() + 
  theme_minimal() + 
  labs(title = 'Priority Selection Status of Sites') 

#' Maximize Dispersion Site Selection
#'
#' Select a subset of sites that maximize spatial dispersion of sites using k-medioids clustering. 
#'
#' @description This function operates on individual points - representing populations, rather than drawing convex hulls or polygons around them to emulate a species range.
#' It is designed for rare species, where individual populations are relatively scarce, e.g. < 100, and have decent location data.
#' It will perform bootstrap re-sampling to better estimate the true range of the extent species, as well as coordinate jittering to better address geo-location quality.
#' After running n of these simulations it will identify the individual networks of sites (co-location) which is the most resilient to these perturbations, and should be less affected by data quality issues.
#' A particular point of this function, relative to the grid based approaches in the package, is that it treats populations as individuals, and allows curators to focus more on  'edges' of species ranges.
#'
#' As arguments it takes the known locations of populations, and will solve for n *priority* collection sites.
#' Along this process it will also generate a priority ranking of all sites, indicating a naive possible order for prioritizing collections; although opportunity should never discard a site.
#' A required input parameter is a column indicating whether a site is a *required*.
#' Required sites (1 - as many as < n_sites) will serve as fixed parameters in the optimization scenario which greatly speed up run time.
#' They can represent: existing collections, collections with a very strong chance of happenging due to a funding agency mechanism, or otherwise a single population closet to the geographic center of the species.
#' Notably the solve will be 'around' this site, hence the solves are not purely theoretical, but linked to a pragmatic element.
#'
#' Theoretically one can substitute a *geographic* distance matrix for an *environmental* distance matrix. See vignette for example. 
#'
#' @param input_data A list with two elements: 'distances' (distance matrix) and 'sites' (data frame of site metadata).
#' @param n_sites The number of sites which you want to select for priority collection.
#' Note that the results will return a rank of prioritization for all sites in the data.
#' @param n_bootstrap Number of bootstrap replicates to perform.
#' @param dropout_prob Probability of dropping non-seed sites in each bootstrap replicate, give how few sites there are generally keep under 0.2. 
#' @param n_local_search_iter Number of local search iterations per restart.
#' @param n_restarts Number of random restarts per bootstrap replicate.
#' @param verbose Whether to print progress information. Will print a message on run settings, and a progress bar for the bootstraps.
#' @examples \dontrun{
#'  library(ggplot2)
#'
#'  n_sites <- 30 # number of known populations
#'  df <- data.frame(
#'    site_id = seq_len(n_sites),
#'    lat = runif(n_sites, 25, 30), # play with these to see elongated results. 
#'    lon = runif(n_sites, -125, -120),
#'    required = FALSE,
#'    coord_uncertainty = 0
#'  )
#'
#' #function relies on at least one required point. here arbitrarily place near geographic center
#'  dists2c <- greatCircleDistance(
#'    median(df$lat), 
#'    median(df$lon), 
#'    df$lat, 
#'    df$lon
#'  )
#'  df[order(dists2c)[1],'required'] <- TRUE
#'  
#'  ## we will simulate coordinate uncertainty on a number of sites.  
#'  uncertain_sites <- sample(setdiff(seq_len(n_sites), which(df$required)), size = min(6, n_sites-3))
#'  df$coord_uncertainty[uncertain_sites] <- runif(length(uncertain_sites), 5000, 100000) # meters
#'  
#'  # the function can take up to take matrices. the first (required) is a geographic distance
#'  # matrix. calculate this with the `greatCircleDistance` fn from the package for consistency. 
#'  # (it will be recalculated during simulations). `sf` gives results in slightly diff units. 
#'  dist_mat <- sapply(1:nrow(df), function(i) {
#'    greatCircleDistance(
#'      df$lat[i], df$lon[i],
#'      df$lat, df$lon
#'    )
#'  })
#'
#'  # the input data is a list, the distance matrix, and the df of actual point locations. 
#'  head(df)
#' 
#'  test_data <- list(distances = dist_mat, sites = df)
#'  rm(dist_mat, df, n_sites, uncertain_sites, dists2c)
#'
#'  # small quick run 
#'    system.time( 
#'      res <- maximizeDispersion(  ## reduce some parameters for faster run. 
#'        input_data = test_data,
#'        n_bootstrap = 500,
#'        n_local_search_iter = 50,
#'        n_restarts = 2
#'      )
#'    )
#'
#' ### first selected 
#'  ggplot(data = res$input_data, 
#'    aes(
#'      x = lon, 
#'      y = lat, 
#'      shape = required, 
#'      size = cooccur_strength,
#'      color = selected
#'      )
#'    ) +
#'    geom_point() + 
#'  #  ggrepel::geom_label_repel(aes(label = site_id), size = 4) + 
#'    theme_minimal() + 
#'    labs(main = 'Priority Selection Status of Sites') 
#'
#' ### order of sampling priority ranking plot.
#'  ggplot(data = res$input_data, 
#'    aes(
#'      x = lon, 
#'      y = lat, 
#'      shape = required, 
#'      size = -sample_rank,
#'      color = sample_rank
#'      )
#'    ) +
#'    geom_point() + 
#' #   ggrepel::geom_label_repel(aes(label = sample_rank), size = 4) +
#'    theme_minimal()   
#' }
#' @export
maximizeDispersion <- function(
    input_data,
    n_sites = 5,
    n_bootstrap = 999,
    dropout_prob = 0.1,
    n_local_search_iter = 100,
    n_restarts = 3,
    verbose = TRUE
) {
  
  distances <- input_data$distances
  sites_df <- input_data$sites
  n_total <- nrow(sites_df)
  
  # Setup required sites and uncertainty
  if (!"required" %in% names(sites_df)) {
    sites_df$required <- FALSE
  }
  seeds <- which(sites_df$required)

  ## identify sites which can be jittered for coord uncertainty. 
  uncertain_idx <- which(   ### minimum distance to initiate this is 10 km
    !is.na(sites_df$coord_uncertainty) & sites_df$coord_uncertainty > 10000
  )
  
  # Initialize tracking matrices
  # we will populate this with the results of each bs iteration. 
  cooccur <- matrix(0L, n_total, n_total)
  
  if (verbose) {
    cat(
      sprintf(
        "Sites: %d | Seeds: %d | Requested: %d | Coord. Uncertain: %d | BS Replicates: %d\n",
        n_total,
        length(seeds),
        n_sites,
        length(uncertain_idx),
        n_bootstrap
      )
    )
  }
  
  # Setup resample parameters
  if (dropout_prob > 0) {
    droppable <- which(!sites_df$required)
    n_drop <- floor(length(droppable) * dropout_prob)
  }
  
  all_solutions <- vector("list", n_bootstrap)
  # Bootstrap loop
  pb <- utils::txtProgressBar(min = 0, max = n_bootstrap, style = 3)
  
  for (b in seq_len(n_bootstrap)) {
    
    # sub-sample the data for the bootstrap. 
    dropped <- sample(droppable, n_drop)
    available_sites <- setdiff(seq_len(n_total), dropped)
    
    # Run single bootstrap iteration
    result <- run_bootstrap_iteration(
      distances = distances,
      sites_df = sites_df,
      n_sites = n_sites,
      seeds = seeds,
      available_sites = available_sites,
      uncertain_idx = uncertain_idx,
      n_local_search_iter = n_local_search_iter,
      n_restarts = n_restarts
    )
    
    # Update co-occurrence matrix
    if (!is.null(result$solution)) {
      idx <- result$solution
      cooccur[idx, idx] <- cooccur[idx, idx] + 1
      
      all_solutions[[b]] <- list(
        sites = sort(result$solution),
        objective = result$objective
      )

    }
    
    utils::setTxtProgressBar(pb, b)
  }
  close(pb)
  ##### identify the set of sites found to be the most stable #####

  # Identify most stable combination
  diag(cooccur) <- 0
  cooccurrence_strength <- rowSums(cooccur)
  
  stability <- data.frame(
    site_id = sites_df$site_id,
    cooccur_strength = cooccurrence_strength,
    is_seed = sites_df$required
  )
  stability <- stability[order(-stability$cooccur_strength), ]
  
  # Find the solution which was returned most often. 
  sol_strings <- sapply(all_solutions, function(x) {
      paste(x$sites, collapse = "-")
  })
  tab <- table(sol_strings)
  best_combo_key <- names(which.max(tab))
  most_stable_frequency <- max(tab) / n_bootstrap
  most_stable_solution <- as.integer(strsplit(best_combo_key, "-")[[1]])
  
  ##########     Prepare output    ###########3

  ## summarize the main takeaways - sample these sites first back onto input data. 
  input_appended <- merge(sites_df, stability, by = 'site_id', all.x = TRUE)
  input_appended <- merge(
    input_appended,
    data.frame(
      site_id = most_stable_solution,
      selected = TRUE
    ),
    by = 'site_id', all.x = TRUE
  )

  # these sites were unselected, convert from NA to bool FALSE. 
  input_appended$selected <- replace(
    input_appended$selected,
    is.na(input_appended$selected),
    FALSE
  )
  
  # rank sites by how often they were selected across all simulations
  # these are the most robust to taxonomic changes, local exinctions, bad ids, 
  # bad data etc. 
  input_appended <- input_appended[
    order(input_appended$cooccur_strength, decreasing = TRUE),
  ]

  # give tied ranks to each co-coccur strengths, support notion that
  # all populations should be sampled, but a priority exists. 
  input_appended$sample_rank <- match(
    -stability$cooccur_strength,
    sort(unique(-stability$cooccur_strength))
  )
  
  list(
    input_data = input_appended,
    selected_sites = most_stable_solution,
    stability_score = most_stable_frequency,
    stability = stability,
    settings = data.frame(
      n_sites = n_sites,
      n_bootstrap = n_bootstrap,
      dropout_prob = dropout_prob,
      n_uncertain = length(uncertain_idx)
    )
  )
}


#' Haversine Distance Calculation
#' 
#' calculate distances between sites (Haversine formula)
#' 
#' @description
#' Calculate geographic distances on the geoid.
#' This results in more accurate distance calculations than a planar system. 
#' Function is mostly used internally by `maximize_dispersion`
#' @param lat1 Double. column holding coords of 'focal' population
#' @param lon1 Double. column holding coords of 'focal' population
#' @param lat2 Double. column holding coords of 'non-focal' population
#' @param lon2 Double. column holding coords of 'non-focal' population
#' 
#' @examples
#' n_sites <- 5 # number of known populations
#'  df <- data.frame(
#'    site_id = seq_len(n_sites),
#'    lat = runif(n_sites, 25, 30), 
#'    lon = runif(n_sites, -125, -120)
#'  )
#' 
#' dist_mat <- sapply(1:nrow(df), function(i) {
#'    greatCircleDistance(
#'      df$lat[i], df$lon[i],
#'      df$lat, df$lon
#'    )
#'  })
#' #head(dist_mat)
#' @export
greatCircleDistance <- function(lat1, lon1, lat2, lon2) {
  # Haversine in km
  to_rad <- pi / 180
  lat1 <- lat1 * to_rad
  lon1 <- lon1 * to_rad
  lat2 <- lat2 * to_rad
  lon2 <- lon2 * to_rad
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  c <- 2 * asin(pmin(1, sqrt(a)))
  r <- 6371 # radius of earth.
  return(r * c)
}



# Single bootstrap iteration
run_bootstrap_iteration <- function(
    distances,
    sites_df,
    n_sites,
    seeds,
    available_sites,
    uncertain_idx,
    n_local_search_iter,
    n_restarts
) {
  
  # Jitter distances if there are uncertain coordinates
  if (length(uncertain_idx) > 0) {
    dist_boot <- update_distances_jitter(
      distances,
      sites_df,
      uncertain_idx
    )
  } else {
    dist_boot <- distances
  }
  
# Single bootstrap iteration
  run_bootstrap_iteration(
    distances,
    sites_df,
    n_sites,
    seeds,
    available_sites,
    uncertain_idx,
    n_local_search_iter,
    n_restarts
) {
  
  # Jitter distances if there are uncertain coordinates
  if (length(uncertain_idx) > 0) {
    dist_boot <- update_distances_jitter(
      distances,
      sites_df,
      uncertain_idx
    )
  } else {
    dist_boot <- distances
  }
  
  best_solution <- NULL
  best_objective <- -Inf
  
  # Multi-restart optimization
  for (restart in seq_len(n_restarts)) {

    # filter sites for opportunities.
    available_idx <- available_sites
    dist_boot_filtered <- dist_boot[available_idx, available_idx, drop = FALSE]

    seeds_in_filtered <- which(available_idx %in% seeds)

    pam_result <- pam_fixed(
      dist_matrix = dist_boot_filtered,
      k = n_sites,
      fixed_ids = seeds_in_filtered
    )
    current_solution <- available_idx[pam_result$medoids]
    current_objective <- -pam_result$total_cost 
    
    # Track best across restarts
    if (current_objective > best_objective) {
      best_solution <- current_solution
      best_objective <- current_objective
    }
  }

  list(
    solution = best_solution,
    objective = best_objective
  )
}
  best_solution <- NULL
  best_objective <- -Inf
  
  # Multi-restart optimization
  for (restart in seq_len(n_restarts)) {

    # filter sites for opportunities.
    available_idx <- available_sites
    dist_boot_filtered <- dist_boot[available_idx, available_idx, drop = FALSE]

    seeds_in_filtered <- which(available_idx %in% seeds)

    pam_result <- pam_fixed(
      dist_matrix = dist_boot_filtered,
      k = n_sites,
      fixed_ids = seeds_in_filtered
    )
    current_solution <- available_idx[pam_result$medoids]
    current_objective <- -pam_result$total_cost 
    
    # Track best across restarts
    if (current_objective > best_objective) {
      best_solution <- current_solution
      best_objective <- current_objective
    }
  }

  list(
    solution = best_solution,
    objective = best_objective
  )
}



function()
    distances,
    sites_df,
    n_sites,
    seeds,
    available_sites,
    uncertain_idx,
    n_local_search_iter,
    n_restarts
) {
  
  # Jitter distances if there are uncertain coordinates
  if (length(uncertain_idx) > 0) {
    dist_boot <- update_distances_jitter(
      distances,
      sites_df,
      uncertain_idx
    )
  } else {
    dist_boot <- distances
  }
  
  best_solution <- NULL
  best_objective <- -Inf
  
  # Multi-restart optimization
  for (restart in seq_len(n_restarts)) {

    # filter sites for opportunities.
    available_idx <- available_sites
    dist_boot_filtered <- dist_boot[available_idx, available_idx, drop = FALSE]

    seeds_in_filtered <- which(available_idx %in% seeds)

    pam_result <- pam_fixed(
      dist_matrix = dist_boot_filtered,
      k = n_sites,
      fixed_ids = seeds_in_filtered
    )
    current_solution <- available_idx[pam_result$medoids]
    current_objective <- -pam_result$total_cost 
    
    # Track best across restarts
    if (current_objective > best_objective) {
      best_solution <- current_solution
      best_objective <- current_objective
    }
  }

  list(
    solution = best_solution,
    objective = best_objective
  )
}

#' Jitter coordinates based on uncertainty (in meters)
#' @keywords internal
#' @noRd
jitter_coords <- function(lat, lon, uncertainty_m) {

  # handle missing / zero uncertainty
  uncertainty_m[is.na(uncertainty_m)] <- 0
  
  theta <- runif(length(lat), 0, 2*pi)
  r <- sqrt(runif(length(lat))) * uncertainty_m

  dx <- r * cos(theta)
  dy <- r * sin(theta)

  dlat <- dy / 111320
  dlon <- dx / (111320 * cos(lat * pi/180))

  list(
    jittered_lat = lat + dlat,
    jittered_lon = lon + dlon
  )
}


#' Jitter coordinates of uncertain sites and update distance matrices
#' @keywords internal
#' @noRd
update_distances_jitter <- function(distances, sites_df, uncertain_idx, env_models = NULL, use_model = FALSE) {

  dist_boot <- distances
  coords_boot <- sites_df[c('lat', 'lon')]
  
  if (length(uncertain_idx) > 0) {
    j <- jitter_coords(
      lat   = coords_boot$lat[uncertain_idx],
      lon   = coords_boot$lon[uncertain_idx],
      uncertainty_m = sites_df$coord_uncertainty[uncertain_idx]
    )

    # only replace if jitter returned non-empty output
  if (length(j$lat) > 0) {
      coords_boot$lat[uncertain_idx] <- j$lat
      coords_boot$lon[uncertain_idx] <- j$lon
    }
  }
    # recompute primary distance matrix
    for (i in uncertain_idx) {
      for (j in seq_len(nrow(sites_df))) {
        d <- greatCircleDistance(
          coords_boot$lat[i],
          coords_boot$lon[i],
          coords_boot$lat[j],
          coords_boot$lon[j]
        )
        dist_boot[i, j] <- d
        dist_boot[j, i] <- d
      }
    }
  
  dist_boot
}


# PAM with fixed medoids
pam_fixed <- function(dist_matrix, k, fixed_ids) {
  dist_matrix <- as.matrix(dist_matrix)
  n <- nrow(dist_matrix)
  
  # Initialize with PAM, then ensure fixed sites are included
  init_pam <- cluster::pam(as.dist(dist_matrix), k)$medoids
  init <- setdiff(as.integer(init_pam), fixed_ids)
  init <- head(init, k - length(fixed_ids))
  
  medoids <- unique(c(fixed_ids, init))
  
  # Cost function: sum of minimum distances to nearest medoid
  get_cost <- function(meds) {
    sum(apply(dist_matrix[, meds, drop = FALSE], 1, min))
  }
  
  # Swap optimization loop
  improved <- TRUE
  while (improved) {
    improved <- FALSE
    for (m in medoids) {
      if (m %in% fixed_ids) next  # skip fixed medoids
      
      for (cand in setdiff(1:n, medoids)) {
        trial <- unique(c(fixed_ids, setdiff(medoids, m), cand))
        
        if (length(trial) != length(medoids)) next
        
        if (get_cost(trial) < get_cost(medoids)) {
          medoids <- trial
          improved <- TRUE
        }
      }
    }
  }
  
  # Final cluster assignments
  clusters <- apply(dist_matrix[, medoids, drop = FALSE], 1, which.min)
  
  list(
    medoids = medoids,
    cluster = clusters,
    total_cost = get_cost(medoids)
  )
}