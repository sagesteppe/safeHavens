#' Ensure that both geographic and climatic distances are in 3D array format
#' @keywords internal
#' @noRd
coerce_to_3d_array <- function(distances, n_sites) {
  if (is.null(distances)) {
    stop("'distances' is NULL")
  }
  # vector
  if (
    is.atomic(distances) &&
    is.vector(distances) &&
    !is.list(distances) &&
    !is.matrix(distances) &&
    !is.array(distances)
  ) {
    if (length(distances) == n_sites * n_sites) {
      mat <- matrix(as.numeric(distances), nrow = n_sites, ncol = n_sites)
      return(array(mat, dim = c(n_sites, n_sites, 1)))
    } else {
      stop(sprintf(
        "distance vector length %d doesn't match %d sites (expected %d)",
        length(distances),
        n_sites,
        n_sites * n_sites
      ))
    }
  }
  # matrix
  if (is.matrix(distances)) {
    if (nrow(distances) != n_sites || ncol(distances) != n_sites) {
      stop("distance matrix dims != n_sites")
    }
    return(array(as.numeric(distances), dim = c(n_sites, n_sites, 1)))
  }
  # array
  if (is.array(distances)) {
    dims <- dim(distances)
    if (length(dims) == 2) {
      if (dims[1] != n_sites || dims[2] != n_sites) {
        stop("2D array dims mismatch")
      }
      return(array(as.numeric(distances), dim = c(n_sites, n_sites, 1)))
    } else if (length(dims) == 3) {
      if (dims[1] != n_sites || dims[2] != n_sites) {
        stop("3D array dims mismatch")
      }
      return(distances)
    } else {
      stop("distances array must be 2D or 3D")
    }
  }
  stop("Unsupported 'distances' type")
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
  n <- nrow(sites_df)
  dist_arr <- coerce_to_3d_array(distances, n)
  K <- dim(dist_arr)[3]
  dist_boot <- array(dist_arr, dim = c(n, n, K))
  coords_boot <- data.frame(
    lat = sites_df$lat,
    lon = sites_df$lon
  )
  
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
      for (j in seq_len(n)) {
        d <- greatCircleDistance(
          coords_boot$lat[i],
          coords_boot$lon[i],
          coords_boot$lat[j],
          coords_boot$lon[j]
        )
        dist_boot[i, j, 1] <- d
        dist_boot[j, i, 1] <- d
      }
    }
  
      # update other layers conditionally
  if (K >= 2 && length(uncertain_idx) > 0) {
    for (k in 2:K) {
      if (use_model && !is.null(env_models) && !is.null(env_models[[k]])) {
        fit <- env_models[[k]]
        for (i in uncertain_idx) {
          for (j in seq_len(n)) {
            new_geo <- dist_boot[i, j, 1]
            # predict via the linear model (conditional expectation)
            pred_env <- stats::predict(fit, newdata = data.frame(geo_vec = new_geo))
            # make sure prediction is finite
            if (!is.finite(pred_env)) {
              pred_env <- dist_boot[i, j, k]  # fallback to original
            }
            dist_boot[i, j, k] <- pred_env
            dist_boot[j, i, k] <- pred_env
          }
        }
      } else {
        # small-n fallback: small multiplicative jitter (~2.5%)
        for (i in uncertain_idx) {
          jitter_factor <- stats::rnorm(1, mean = 1, sd = 0.025)
          dist_boot[i, , k] <- dist_boot[i, , k] * jitter_factor
          dist_boot[, i, k] <- dist_boot[, i, k] * jitter_factor
        }
      }
    }
  }
  dist_boot
}

# Greedy initialization
#' @keywords internal
#' @noRd
greedy_initialize <- function(
    dist_matrix,
    n_sites,
    seeds = integer(0),
    objective = "sum"
) {
  n_total <- nrow(dist_matrix)
  selected <- as.integer(seeds)
  if (length(selected) > n_sites) {
    selected <- selected[1:n_sites]
  }
  candidates <- setdiff(seq_len(n_total), selected)
  
  while (length(selected) < n_sites && length(candidates) > 0) {
    best_site <- NULL
    best_obj <- -Inf
    for (site in candidates) {
      test_selected <- c(selected, site)
      if (objective == "sum") {
        obj <- calc_objective_sum(dist_matrix, as.integer(test_selected))
      } else {
        obj <- calc_objective_maxmin(dist_matrix, as.integer(test_selected))
      }
      if (obj > best_obj) {
        best_obj <- obj
        best_site <- site
      }
    }
    if (is.null(best_site)) {
      break
    }
    selected <- c(selected, best_site)
    candidates <- setdiff(candidates, best_site)
  }
  as.integer(selected)
}

# Variance-penalized sum objective
#' @keywords internal
#' @noRd
calc_objective_sum_var <- function(dist_mat, selected, lambda_var) {
  n <- length(selected)
  if (n < 2) {
    return(0)
  }
  
  dists <- dist_mat[selected, selected]
  dists <- dists[lower.tri(dists)]
  
  obj_disp <- sum(dists)
  obj_var <- if (length(dists) >= 2) -var(dists) else 0
  
  return(obj_disp + lambda_var * obj_var)
}

# Wrapper for greedy initialization using new objective
#' @keywords internal
#' @noRd
greedy_initialize_var <- function(dist_matrix, n_sites, seeds, lambda_var) {
  n_total <- nrow(dist_matrix)
  selected <- seeds
  candidates <- setdiff(1:n_total, seeds)
  
  while (length(selected) < n_sites && length(candidates) > 0) {
    best_site <- NULL
    best_obj <- -Inf
    
    for (site in candidates) {
      test_selected <- c(selected, site)
      obj <- calc_objective_sum_var(dist_matrix, test_selected, lambda_var)
      if (obj > best_obj) {
        best_obj <- obj
        best_site <- site
      }
    }
    
    selected <- c(selected, best_site)
    candidates <- setdiff(candidates, best_site)
  }
  
  return(as.integer(selected))
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

#' Maximize Dispersion Site Selection
#'
#' Select a subset of sites that maximize spatial dispersion using a bootstrapped local search algorithm, combining both maximizing distance and minimizing variance.
#'
#' @description This function operates on individual points, rather than drawing convex hulls or polygons around them.
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
#'
#' @param input_data A list with two elements: 'distances' (distance matrix or array) and 'sites' (data frame of site metadata).
#' @param lambda_var Essentially a smoothing parameter that controls the trade-off between maximizing dispersion and minimizing variance in pairwise distances among selected sites.
#' higher values prioritize variance reduction more strongly.
#' We recommend checking stops between 0.05 and 0.3 to see what works best for your data.
#' Also check up at 0.5+, when getting started, to get a feel for how strong the variance reduction penalization can become. 
#' @param n_sites The number of sites which you want to select for priority collection.
#' Note that the results will return a rank of prioritization for all sites in the data.
#' @param weight_1 Weights for combining multiple distance matrices (if provided).
#' weight_1 is for the *geographic distance* matrix
#' @param weight_2 Weights for combining multiple distance matrices (if provided).
#' weight_2 is for the *climatic distance* matrix (if provided).
#' @param n_bootstrap Number of bootstrap replicates to perform.
#' @param dropout_prob Probability of dropping non-seed sites in each bootstrap replicate.
#' @param objective Objective function to optimize: "sum" (dispersion sum with variance penalty) or "maxmin" (maximize minimum distance).
#' @param n_local_search_iter Number of local search iterations per restart.
#' @param n_restarts Number of random restarts per bootstrap replicate.
#' @param seed Random seed for reproducibility.
#' @param verbose Whether to print progress information. Will print a message on run settings, and a progress bar for the bootstraps.
#' @examples
#' 
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
#'  df$coord_uncertainty[uncertain_sites] <- runif(length(uncertain_sites), 1000, 10000) # meters
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
#'  # small quick run (fast); timing for R package / CRAN status. 
#'    system.time( 
#'      res <- maximizeDispersion(  ## reduce some parameters for faster run. 
#'        input_data = test_data,
#'        lambda_var = 0.2,
#'        n_bootstrap = 500,
#'        objective = "maxmin",
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
#' 
#' @export
maximizeDispersion <- function(
    input_data,
    lambda_var = 0.15,
    n_sites = 5,
    weight_1 = 1.0,
    weight_2 = 0.0,
    n_bootstrap = 999,
    dropout_prob = 0.15,
    objective = c("sum", "maxmin"),
    n_local_search_iter = 100,
    n_restarts = 3,
    seed = NULL,
    verbose = TRUE
) {
  objective <- match.arg(objective)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  distances <- input_data$distances
  sites_df <- input_data$sites
  n_total <- nrow(sites_df)
  cooccur <- matrix(0, n_total, n_total)
  
  if (n_total <= 0) {
    stop("No sites provided")
  }
  
  # coerce to 3D
  distances <- coerce_to_3d_array(distances, n_total)
  K <- dim(distances)[3]

  
  # combine first two matrices by weights
  dist_combined <- distances[,, 1] * weight_1
  if (dim(distances)[3] >= 2 && weight_2 > 0) {
    dist_combined <- dist_combined + distances[,, 2] * weight_2
  }
  
  # set up matrix so that only non-required sites can be dropped from the permutations.
  if (!"required" %in% names(sites_df)) {
    sites_df$required <- FALSE
  }
  sites_df$required <- as.logical(sites_df$required)
  seeds <- which(sites_df$required)
  uncertain_idx <- which(
    !is.na(sites_df$coord_uncertainty) & sites_df$coord_uncertainty > 0
  )
  
  #####################################################################
  # FIT env_models this is used to jitter the environmental distances between sites
  #####################################################################
  use_model <- (K >= 2) && (n_total > 25)
  env_models <- NULL
  if (use_model) {
    geo_vec <- distances[,,1][lower.tri(distances[,,1])]
    env_models <- vector("list", K)
    env_models[[1]] <- NA
    for (k in 2:K) {
      env_vec <- distances[,,k][lower.tri(distances[,,k])]
      # safe lm fit: if env_vec or geo_vec invalid, set NA
      if (length(env_vec) == length(geo_vec) && length(geo_vec) > 1 && all(is.finite(geo_vec)) && all(is.finite(env_vec))) {
        env_models[[k]] <- tryCatch(lm(env_vec ~ geo_vec), error = function(e) NULL)
      } else {
        env_models[[k]] <- NULL
      }
    }
  }

  cooccur <- matrix(0L, n_total, n_total)
  selection_counts <- integer(n_total)
  all_solutions <- list()
  solution_counter <- 1

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
  
  ##############################################
  ###                bootstrap              ###
  ##############################################
  if (dropout_prob > 0) {
    droppable <- setdiff(seq_len(n_total), seeds)
    n_drop <- floor(length(droppable) * dropout_prob)
    should_dropout <- (n_drop > 0 && length(droppable) >= n_drop)
  }
  
  pb <- utils::txtProgressBar(min = 0, max = n_bootstrap, style = 3)
  
  for (b in seq_len(n_bootstrap)) {
    ## bs begins.
    available_sites <- seq_len(n_total)
    
    # dropout (non-seed only)
    if (should_dropout) {
      dropped <- sample(droppable, n_drop)
      available_sites <- setdiff(available_sites, dropped)
    }
    
    # jitter distances if needed
    if (length(uncertain_idx) > 0) {
      dist_boot_array <- update_distances_jitter(
        distances,
        sites_df,
        uncertain_idx, 
        env_models, 
        use_model
      )
      dist_boot <- dist_boot_array[,, 1] * weight_1
      if (K >= 2 && weight_2 > 0) {
        dist_boot <- dist_boot + dist_boot_array[,,2] * weight_2
      }
    } else {
      dist_boot <- dist_combined
    }
    
    best_solution <- NULL
    best_objective <- -Inf
    
    for (restart in seq_len(n_restarts)) {
      current_solution <- greedy_initialize_var(
        dist_boot,
        n_sites,
        seeds,
        lambda_var = lambda_var
      )
      # ensure selection length
      if (length(current_solution) < n_sites) {
        # try to pad with available_sites
        extra <- setdiff(available_sites, current_solution)
        if (length(extra) > 0) {
          needed <- n_sites - length(current_solution)
          current_solution <- c(current_solution, utils::head(extra, needed))
        }
      }
      
      swappable_solution <- setdiff(current_solution, seeds)
      candidates <- setdiff(available_sites, current_solution)
      
      if (length(candidates) > 0 && length(swappable_solution) > 0) {
        res <- local_search_swap(
          dist_boot,
          as.integer(swappable_solution),
          as.integer(candidates),
          objective,
          n_local_search_iter,
          lambda_var
        )
        sol <- as.integer(res$selected)
        objv <- as.numeric(res$objective)
        if (objv > best_objective) {
          best_objective <- objv
          best_solution <- c(seeds, as.integer(res$selected))
        }
      } else {
        # compute objective for current_solution
        if (length(current_solution) > 0) {
          if (objective == "sum") {
            curv <- calc_objective_sum_var(
              dist_boot,
              as.integer(current_solution),
              lambda_var = lambda_var
            )
          } else {
            curv <- calc_objective_maxmin(
              dist_boot,
              as.integer(current_solution)
            )
          }
          if (curv > best_objective) {
            best_objective <- curv
            best_solution <- as.integer(current_solution)
          }
        }
      }
    } # restarts
    
    if (!is.null(best_solution)) {
      idx <- best_solution
      cooccur[idx, idx] <- cooccur[idx, idx] + 1
      
      all_solutions[[solution_counter]] <- list(
        sites = sort(best_solution),
        objective = best_objective
      )
      solution_counter <- solution_counter + 1
    }
    
    utils::setTxtProgressBar(pb, b)
  } # end  bootstrap
  close(pb) # progress bar 
  
  
  ###########################################
  ### Identify the most stable combination ###
  ## the combination that appears most frequently across bootstraps. ##
  ###########################################
  # remove diagonal — marginal selection is irrelevant
  diag(cooccur) <- 0
  
  # co-selection strength = sum of pairwise affinities
  cooccurrence_strength <- rowSums(cooccur)
  
  # boost seeds so they always appear in top ranks
  if (length(seeds) > 0) {
    cooccurrence_strength[seeds] <- max(cooccurrence_strength) + 1
  }
  
  stability <- data.frame(
    site_id = sites_df$site_id,
    cooccur_strength = cooccurrence_strength,
    is_seed = sites_df$required
  )
  
  # rank by co-selection
  stability <- stability[order(-stability$cooccur_strength), ]
  
  if (length(all_solutions) == 0) {
    most_stable_solution <- rep(NA, n_sites)
    most_stable_frequency <- 0
  } else {
    # Convert combos to unique keys
    sol_strings <- sapply(all_solutions, function(x) {
      paste(x$sites, collapse = "-")
    })
    
    # Count frequency of each unique combination
    tab <- table(sol_strings)
    
    # Most frequent combo across bootstraps
    best_combo_key <- names(which.max(tab))
    most_stable_frequency <- max(tab) / n_bootstrap
    most_stable_solution <- as.integer(strsplit(best_combo_key, "-")[[1]])
  }
  
  
  ##############################################
  ###             return objects            ####
  ##############################################
  ## input data with a couple columns added on, probably all most users want.
  input_appended = merge(sites_df, stability, on = 'site_id', how = 'left')
  input_appended = merge(
    input_appended,
    data.frame(
      site_id = most_stable_solution,
      selected = T
    ),
    on = 'site_id',
    how = 'left',
    all.x = T
  )
  input_appended$selected <- replace(
    input_appended$selected,
    is.na(input_appended$selected),
    FALSE
  )  
  
  ### add on the relative selection/collection priority ranks for the populations. 
  input_appended <- input_appended[
    order(input_appended$cooccur_strength, decreasing = TRUE),
  ]
  input_appended$sample_rank <- match(
    -stability$cooccur_strength,
    sort(unique(-stability$cooccur_strength))
  )
  
  # return objects back to the user.
  list(
    input_data = input_appended,
    selected_sites = most_stable_solution,
    stability_score = most_stable_frequency,
    stability = stability,
    settings = data.frame(
      n_sites = n_sites,
      n_bootstrap = n_bootstrap,
      objective = objective,
      lambda = lambda_var,
      dropout_prob = dropout_prob,
      n_uncertain = length(uncertain_idx)
    )
  )
}