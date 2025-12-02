# Example usage
set.seed(20)
library(ggplot2)

n_sites <- 30 # number of known populations
df <- data.frame(
  site_id = seq_len(n_sites),
  lat = runif(n_sites, 25, 30), # play with these to see elongated results. 
  lon = runif(n_sites, -125, -120),
  required = FALSE,
  coord_uncertainty = 0
)

## we will simulate coordinate uncertainty on a number of sites.  
uncertain_sites <- sample(setdiff(seq_len(n_sites), which(df$required)), size = min(6, n_sites-3))
df$coord_uncertainty[uncertain_sites] <- runif(length(uncertain_sites), 1000, 10000) # meters

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
  distances = dist_mat,
  sites = df
  )

resin <- maximizeDispersion(dat, n_sites = 6, algorithm = 'pam', n_bootstrap=99, dropout_prob = 0.2)

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
#  ggrepel::geom_label_repel(aes(label = site_id), size = 4) + 
  theme_minimal() + 
  labs(title = 'Priority Selection Status of Sites') 


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

# Single bootstrap iteration
run_bootstrap_iteration <- function(
    dist_combined,
    distances,
    sites_df,
    n_sites,
    seeds,
    available_sites,
    uncertain_idx,
    env_models,
    use_model,
    weight_1,
    weight_2,
    K,
    lambda_var,
    objective,
    n_local_search_iter,
    n_restarts,
    algorithm = c("greedy", "pam")
) {
  algorithm <- match.arg(algorithm)
  
  # Jitter distances if there are uncertain coordinates
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
      dist_boot <- dist_boot + dist_boot_array[,, 2] * weight_2
    }
  } else {
    dist_boot <- dist_combined
  }
  
  best_solution <- NULL
  best_objective <- -Inf
  
  # Multi-restart optimization
  for (restart in seq_len(n_restarts)) {
    
    if (algorithm == "pam") {

      available_idx <- available_sites
      dist_boot_filtered <- dist_boot[available_idx, available_idx, drop = FALSE]

      seeds_in_filtered <- which(available_idx %in% seeds)

      # Use PAM algorithm
      pam_result <- pam_fixed(
        dist_matrix = dist_boot_filtered,
        k = n_sites,
        fixed_ids = seeds_in_filtered
      )
      current_solution <- available_idx[pam_result$medoids]
      # For PAM, use negative cost as objective (minimization -> maximization)
      objv <- -pam_result$total_cost
      
    } else {
      # Use greedy initialization
      current_solution <- greedy_initialize_var(
        dist_boot,
        n_sites,
        seeds,
        lambda_var = lambda_var
      )
      
      # Ensure solution has correct length
      if (length(current_solution) < n_sites) {
        extra <- setdiff(available_sites, current_solution)
        if (length(extra) > 0) {
          needed <- n_sites - length(current_solution)
          current_solution <- c(current_solution, utils::head(extra, needed))
        }
      }
      
      swappable_solution <- setdiff(current_solution, seeds)
      candidates <- setdiff(available_sites, current_solution)
      
      # Local search with swaps
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
        current_solution <- c(seeds, sol)
      } else {
        # Compute objective for current solution
        if (length(current_solution) > 0) {
          if (objective == "sum") {
            objv <- calc_objective_sum_var(
              dist_boot,
              as.integer(current_solution),
              lambda_var = lambda_var
            )
          } else {
            objv <- calc_objective_maxmin(
              dist_boot,
              as.integer(current_solution)
            )
          }
        } else {
          objv <- -Inf
        }
      }
    }
    
    # Track best solution across restarts
    if (objv > best_objective) {
      best_objective <- objv
      best_solution <- as.integer(current_solution)
    }
  }
  
  list(
    solution = best_solution,
    objective = best_objective
  )
}

# Main function with abstracted bootstrap
maximizeDispersion <- function(
    input_data,
    lambda_var = 0.15,
    n_sites = 5,
    weight_1 = 1.0,
    weight_2 = 0.0,
    n_bootstrap = 999,
    dropout_prob = 0.15,
    objective = c("sum", "maxmin"),
    algorithm = c("greedy", "pam"),
    n_local_search_iter = 100,
    n_restarts = 3,
    seed = NULL,
    verbose = TRUE
) {
  objective <- match.arg(objective)
  algorithm <- match.arg(algorithm)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  distances <- input_data$distances
  sites_df <- input_data$sites
  n_total <- nrow(sites_df)
  print(n_total)
  
  if (n_total <= 0) {
    stop("No sites provided")
  }
  
  # Coerce to 3D array
  distances <- coerce_to_3d_array(distances, n_total)
  K <- dim(distances)[3]
  
  # Combine first two matrices by weights
  dist_combined <- distances[,, 1] * weight_1
  if (K >= 2 && weight_2 > 0) {
    dist_combined <- dist_combined + distances[,, 2] * weight_2
  }
  
  # Setup required sites and uncertainty
  if (!"required" %in% names(sites_df)) {
    sites_df$required <- FALSE
  }
  sites_df$required <- as.logical(sites_df$required)
  seeds <- which(sites_df$required)
  uncertain_idx <- which(
    !is.na(sites_df$coord_uncertainty) & sites_df$coord_uncertainty > 0
  )
  
  # Fit environmental models for jittering
  use_model <- (K >= 2) && (n_total > 25)
  env_models <- NULL
  if (use_model) {
    geo_vec <- distances[,, 1][lower.tri(distances[,, 1])]
    env_models <- vector("list", K)
    env_models[[1]] <- NA
    for (k in 2:K) {
      env_vec <- distances[,, k][lower.tri(distances[,, k])]
      if (length(env_vec) == length(geo_vec) && length(geo_vec) > 1 && 
          all(is.finite(geo_vec)) && all(is.finite(env_vec))) {
        env_models[[k]] <- tryCatch(lm(env_vec ~ geo_vec), error = function(e) NULL)
      } else {
        env_models[[k]] <- NULL
      }
    }
  }
  
  # Initialize tracking matrices
  cooccur <- matrix(0L, n_total, n_total)
  all_solutions <- list()
  
  if (verbose) {
    cat(
      sprintf(
        "Sites: %d | Seeds: %d | Requested: %d | Coord. Uncertain: %d | BS Replicates: %d | Algorithm: %s\n",
        n_total,
        length(seeds),
        n_sites,
        length(uncertain_idx),
        n_bootstrap,
        algorithm
      )
    )
  }
  
  # Setup dropout parameters
  should_dropout <- FALSE
  if (dropout_prob > 0) {
    droppable <- setdiff(seq_len(n_total), seeds)
    n_drop <- floor(length(droppable) * dropout_prob)
    should_dropout <- (n_drop > 0 && length(droppable) >= n_drop)
  }
  
  # Bootstrap loop
  pb <- utils::txtProgressBar(min = 0, max = n_bootstrap, style = 3)
  
  for (b in seq_len(n_bootstrap)) {
    available_sites <- seq_len(n_total)
    
    # Apply dropout
    if (should_dropout) {
      dropped <- sample(droppable, n_drop)
      available_sites <- setdiff(available_sites, dropped)
    }
    
    # Run single bootstrap iteration
    result <- run_bootstrap_iteration(
      dist_combined = dist_combined,
      distances = distances,
      sites_df = sites_df,
      n_sites = n_sites,
      seeds = seeds,
      available_sites = available_sites,
      uncertain_idx = uncertain_idx,
      env_models = env_models,
      use_model = use_model,
      weight_1 = weight_1,
      weight_2 = weight_2,
      K = K,
      lambda_var = lambda_var,
      objective = objective,
      n_local_search_iter = n_local_search_iter,
      n_restarts = n_restarts,
      algorithm = algorithm
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
  
  # Identify most stable combination
  diag(cooccur) <- 0
  cooccurrence_strength <- rowSums(cooccur)
  
  if (length(seeds) > 0) {
    cooccurrence_strength[seeds] <- max(cooccurrence_strength) + 1
  }
  
  stability <- data.frame(
    site_id = sites_df$site_id,
    cooccur_strength = cooccurrence_strength,
    is_seed = sites_df$required
  )
  stability <- stability[order(-stability$cooccur_strength), ]
  
  # Find most frequent solution
  if (length(all_solutions) == 0) {
    most_stable_solution <- rep(NA, n_sites)
    most_stable_frequency <- 0
  } else {
    sol_strings <- sapply(all_solutions, function(x) {
      paste(x$sites, collapse = "-")
    })
    tab <- table(sol_strings)
    best_combo_key <- names(which.max(tab))
    most_stable_frequency <- max(tab) / n_bootstrap
    most_stable_solution <- as.integer(strsplit(best_combo_key, "-")[[1]])
  }
  
  # Prepare output
  input_appended <- merge(sites_df, stability, by = 'site_id', all.x = TRUE)
  input_appended <- merge(
    input_appended,
    data.frame(
      site_id = most_stable_solution,
      selected = TRUE
    ),
    by = 'site_id',
    all.x = TRUE
  )
  input_appended$selected <- replace(
    input_appended$selected,
    is.na(input_appended$selected),
    FALSE
  )
  
  input_appended <- input_appended[
    order(input_appended$cooccur_strength, decreasing = TRUE),
  ]
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
      objective = objective,
      algorithm = algorithm,
      lambda = lambda_var,
      dropout_prob = dropout_prob,
      n_uncertain = length(uncertain_idx)
    )
  )
}
