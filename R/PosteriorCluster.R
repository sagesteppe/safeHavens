#' Cluster environmental space using posterior beta draws
#'
#' @description
#' Implements option D of the Bayesian clustering workflow: draws K sets of
#' beta coefficients from the posterior, produces a weighted environmental
#' raster for each draw, samples points from each, and records how often pairs
#' of sample points are assigned to the same cluster. The resulting co-occurrence
#' probability matrix is then used to derive a consensus clustering.
#'
#' Cells that co-cluster with probability near 1 across all draws are robustly
#' similar in environmental space. Cells near boundaries that shift across draws
#' — where the model is uncertain which variables matter most — will have
#' intermediate co-occurrence probability. This boundary-instability surface is
#' returned separately as `StabilityRaster` and is directly meaningful for
#' targeting collection effort.
#'
#' @param model A `brmsfit` object from [bayesianSDM()].
#' @param predictors A `SpatRaster` stack (`$Predictors` from [bayesianSDM()]).
#' @param f_rasts The rasters output from [bayesianSDM()], used to define the
#'   spatial mask and domain. Specifically uses `$RasterPredictions`.
#' @param pred_mat Data frame or matrix (`$PredictMatrix` from [bayesianSDM()]).
#' @param training_data `sf` object (`$TrainData` from [bayesianSDM()]).
#' @param n_draws Integer. Number of posterior beta draws to cluster over.
#'   Defaults to `100`. Values between 50–200 give stable co-occurrence
#'   matrices for typical SDM datasets (< 2000 obs, < 20 vars).
#' @param n Integer. Number of clusters per draw. Passed to
#'   [EnvironmentalBasedSample()] internals.
#' @param n_pts Integer. Points to sample per draw for clustering. Defaults to
#'   `500`. The same points are used across all draws (fixed spatial frame)
#'   so co-occurrence is directly comparable.
#' @param lyr Character. Layer name in `f_rasts` to use as the spatial domain
#'   mask. Defaults to `"occurrence_prob_mean"` (the [bayesianSDM()] output
#'   name).
#' @param planar_proj Numeric EPSG or proj4 string for planar projection.
#' @param coord_wt Numeric. Coordinate weight, as in [EnvironmentalBasedSample()].
#'   Defaults to `2.5`.
#' @param consensus_method Character. How to derive the final clustering from
#'   the co-occurrence matrix. One of:
#'   \describe{
#'     \item{`"hierarchical"`}{Average-linkage hierarchical clustering on
#'       `(1 - co_occurrence)` as a dissimilarity. Stable and interpretable.
#'       Default.}
#'     \item{`"pam"`}{Partitioning Around Medoids on the dissimilarity matrix.
#'       More robust to outliers than hierarchical.}
#'   }
#' @param beta_draws Optional matrix of pre-drawn posterior samples (rows =
#'   draws, cols = parameters, as from [brms::as_draws_matrix()]). If `NULL`,
#'   draws are taken internally. Supply this to reuse draws across multiple
#'   calls.
#' @param seed Integer. Random seed. Defaults to `42`.
#'
#' @returns A list:
#' \describe{
#'   \item{`Geometry`}{`sf` polygon layer of consensus clusters (compatible
#'     with [EnvironmentalBasedSample()] `$Geometry` output).}
#'   \item{`CoOccurrenceMatrix`}{Numeric matrix (n_pts × n_pts). Entry [i,j]
#'     is the proportion of draws in which points i and j were assigned to the
#'     same cluster.}
#'   \item{`StabilityRaster`}{`SpatRaster` of per-cell cluster-assignment
#'     stability, predicted via KNN regression from the sample point stability
#'     scores. Values in [0,1]; low values flag environmentally ambiguous
#'     boundary regions.}
#'   \item{`ConsensusRaster`}{`SpatRaster` of final consensus (rank1) cluster
#'     IDs, geographically reordered.}
#'   \item{`Rank2Raster`}{`SpatRaster` of second-most-frequent cluster
#'     assignments across draws. Where a point was always assigned to one
#'     cluster, rank2 = rank1.}
#'   \item{`Rank3Raster`}{`SpatRaster` of third-most-frequent cluster
#'     assignments. Where fewer than 3 unique clusters were assigned, rank3 =
#'     rank1.}
#'   \item{`Top3Lookup`}{Data frame with columns: `point_id`, `x`, `y`,
#'     `rank1_cluster`, `rank1_pct`, `rank2_cluster`, `rank2_pct`,
#'     `rank3_cluster`, `rank3_pct`, `uncertainty` (= 100 - rank1_pct).
#'     Tabular summary of cluster assignment frequencies at sample points.}
#'   \item{`DrawClusterings`}{Integer matrix (n_pts × n_draws) of per-draw
#'     cluster assignments. Retained for diagnostics.}
#'   \item{`SamplePoints`}{`sf` object of the fixed sample points used.}
#'   \item{`KNN_Cluster`}{Trained KNN model for rank1 cluster assignment. Can
#'     be reused to classify new points.}
#'   \item{`KNN_Rank2`}{Trained KNN model for rank2 clusters.}
#'   \item{`KNN_Rank3`}{Trained KNN model for rank3 clusters.}
#'   \item{`KNN_Stability`}{Trained KNN regression model for stability scores.
#'     Can be reused to predict stability at new locations.}
#' }
#'
#' @seealso [RescaleRasters_bayes()], [EnvironmentalBasedSample()],
#'   [bayesianSDM()]
#' @export
PosteriorCluster <- function(
    model,
    predictors,
    f_rasts,
    pred_mat,
    training_data,
    n_draws          = 100,
    n                = 10,
    n_pts            = 500,
    lyr              = "occurrence_prob_mean",
    planar_proj,
    coord_wt         = 2.5,
    consensus_method = c("hierarchical", "pam"),
    beta_draws       = NULL,
    seed             = 42
) {
  set.seed(seed)
  consensus_method <- match.arg(consensus_method)

  # ── 1. Identify environmental predictor names ────────────────────────────────
  # Drop GP coordinate columns; keep only named fixed effects
  fe_names  <- rownames(brms::fixef(model))
  fe_names  <- sub("^b_", "", fe_names)
  fe_names  <- fe_names[!grepl("^(Intercept|sgp|sdgp|lscale)", fe_names)]
  env_vars  <- fe_names[fe_names %in% names(predictors)]

  if (length(env_vars) == 0) {
    stop("No environmental fixed effects matched raster layer names.")
  }

  # ── 2. Draw posterior beta samples ────────────────────────────────────────────
  if (is.null(beta_draws)) {
    message(sprintf("Drawing %d posterior beta samples ...", n_draws))
    beta_draws <- extract_beta_draws(model, env_vars, n_draws)
  }
  # beta_draws: n_draws × length(env_vars), named columns

  # ── 3. Compute training-data moments (fixed across draws) ────────────────────
  sdN <- function(x) sqrt(mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
  pm  <- as.data.frame(pred_mat)[, env_vars, drop = FALSE]
  var_mu <- vapply(env_vars, function(v) mean(pm[[v]], na.rm = TRUE), numeric(1))
  var_sd <- vapply(env_vars, function(v) sdN(pm[[v]]),                 numeric(1))
  zero_sd <- var_sd == 0
  if (any(zero_sd)) {
    warning(sprintf(
      "Variable(s) with zero SD excluded from draw clustering: %s",
      paste(env_vars[zero_sd], collapse = ", ")
    ))
    env_vars  <- env_vars[!zero_sd]
    var_mu    <- var_mu[!zero_sd]
    var_sd    <- var_sd[!zero_sd]
    beta_draws <- beta_draws[, env_vars, drop = FALSE]
  }

  # ── 4. Fix sample points across all draws ────────────────────────────────────
  # Points are sampled once from the SDM prediction surface and held constant.
  # This is essential: co-occurrence is only meaningful when comparing the
  # same spatial locations across draws.
  message("Sampling fixed point set from prediction surface ...")
  mask_rast <- f_rasts[[lyr]]
  sample_pts <- terra::spatSample(
    mask_rast,
    size      = n_pts,
    method    = "random",
    as.points = TRUE,
    na.rm     = TRUE
  )

  # Extract predictor values at these fixed points
  pt_env <- terra::extract(predictors[[env_vars]], sample_pts, ID = FALSE)
  complete_rows <- stats::complete.cases(pt_env)
  sample_pts    <- sample_pts[complete_rows, ]
  pt_env        <- pt_env[complete_rows, , drop = FALSE]
  n_pts_actual  <- nrow(pt_env)

  message(sprintf(
    "  %d of %d requested points have complete predictor data.",
    n_pts_actual, n_pts
  ))

  # ── 5. Add weighted coordinates to point matrix ───────────────────────────────
  # Coordinates are added once (they don't vary across beta draws)
  pt_env_coords <- add_coord_weights_to_points(
    sample_pts, pt_env, coord_wt, env_vars
  )

  # ── 6. Cluster over posterior draws ──────────────────────────────────────────
  message(sprintf("Clustering over %d posterior draws ...", n_draws))
  draw_clusterings <- matrix(
    NA_integer_,
    nrow = n_pts_actual,
    ncol = n_draws
  )

  for (d in seq_len(n_draws)) {
    if (d %% 25 == 0) message(sprintf("  Draw %d / %d", d, n_draws))

    betas_d <- beta_draws[d, ]  # named vector, length = env_vars

    # Rescale point matrix by this draw's betas
    pt_scaled <- rescale_points_by_betas(pt_env, env_vars, var_mu, var_sd, betas_d)

    # Append (already-weighted) coordinates — same across draws
    pt_full <- cbind(
      pt_scaled,
      pt_env_coords[, c("coord_x_w", "coord_y_w"), drop = FALSE]
    )
    pt_full <- pt_full[stats::complete.cases(pt_full), , drop = FALSE]

    # Hierarchical clustering on Euclidean distances
    dist_d  <- stats::dist(pt_full, method = "euclidean")
    hc_d    <- stats::hclust(dist_d, method = "ward.D2")
    draw_clusterings[seq_len(nrow(pt_full)), d] <- stats::cutree(hc_d, k = n)
  }

  # ── 7. Build co-occurrence matrix ────────────────────────────────────────────
  message("Computing co-occurrence matrix ...")
  co_mat <- build_cooccurrence_matrix(draw_clusterings, n_pts_actual, n_draws)

  # ── 8. Consensus clustering ───────────────────────────────────────────────────
  message("Deriving consensus clustering ...")
  diss_mat         <- stats::as.dist(1 - co_mat)
  consensus_labels <- switch(
    consensus_method,
    hierarchical = {
      hc <- stats::hclust(diss_mat, method = "average")
      stats::cutree(hc, k = n)
    },
    pam = {
      cluster::pam(diss_mat, k = n, diss = TRUE)$clustering
    }
  )

  # ── 9. Stability surface ──────────────────────────────────────────────────────
  # For each point, stability = mean co-occurrence with all other points
  # sharing its consensus cluster label (its "within-cluster cohesion")
  stability_scores <- compute_stability_scores(co_mat, consensus_labels)
  
  # ── 9b. Top-3 cluster rankings ───────────────────────────────────────────────
  message("Computing top-3 cluster assignments per point ...")
  top3_results <- compute_top3_clusters(draw_clusterings, n_pts_actual, n_draws)

  # ── 10. Project clusters back to raster ──────────────────────────────────────
  message("Projecting consensus clusters to raster ...")
  
  # Compute posterior mean betas
  mean_betas <- colMeans(beta_draws)
  
  rast_list <- project_consensus_to_raster(
    sample_pts       = sample_pts,
    consensus_labels = consensus_labels,
    stability_scores = stability_scores,
    top3_labels      = top3_results$labels_matrix,  # n_pts × 3 matrix
    predictors       = predictors,
    env_vars         = env_vars,
    var_mu           = var_mu,
    var_sd           = var_sd,
    mean_betas       = mean_betas,
    coord_wt         = coord_wt,
    mask_rast        = mask_rast,
    planar_proj      = planar_proj
  )

  # ── 11. Geographic reordering ─────────────────────────────────────────────────
  reordered     <- reorder_clusters_geographically(rast_list$cluster_raster)
  final_raster  <- terra::project(reordered$raster, terra::crs(f_rasts))
  cluster_vects <- sf::st_transform(reordered$vectors, terra::crs(f_rasts))
  
  # Project stability raster back to original CRS
  stability_final <- terra::project(rast_list$stability_raster, terra::crs(f_rasts))
  
  # Project top-3 rasters
  rank2_final <- terra::project(rast_list$rank2_raster, terra::crs(f_rasts))
  rank3_final <- terra::project(rast_list$rank3_raster, terra::crs(f_rasts))
  
  # Build tabular lookup with sample point coordinates and top-3 assignments
  sample_pts_sf <- sf::st_as_sf(sample_pts)
  coords_df <- as.data.frame(sf::st_coordinates(sample_pts_sf))
  
  top3_lookup <- data.frame(
    point_id      = seq_len(nrow(top3_results$labels_matrix)),
    x             = coords_df$X,
    y             = coords_df$Y,
    rank1_cluster = top3_results$labels_matrix[, 1],
    rank1_pct     = round(top3_results$pct_matrix[, 1], 1),
    rank2_cluster = top3_results$labels_matrix[, 2],
    rank2_pct     = round(top3_results$pct_matrix[, 2], 1),
    rank3_cluster = top3_results$labels_matrix[, 3],
    rank3_pct     = round(top3_results$pct_matrix[, 3], 1),
    uncertainty   = round(100 - top3_results$pct_matrix[, 1], 1)  # 1 - rank1_pct
  )

  list(
    Geometry           = cluster_vects,
    CoOccurrenceMatrix = co_mat,
    StabilityRaster    = stability_final,
    ConsensusRaster    = final_raster,
    Rank2Raster        = rank2_final,
    Rank3Raster        = rank3_final,
    Top3Lookup         = top3_lookup,
    DrawClusterings    = draw_clusterings,
    SamplePoints       = sample_pts_sf,
    KNN_Cluster        = rast_list$knn_cluster,
    KNN_Rank2          = rast_list$knn_rank2,
    KNN_Rank3          = rast_list$knn_rank3,
    KNN_Stability      = rast_list$knn_stability
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

#' Train a regression KNN for predicting continuous stability scores
#'
#' Uses caret to train a k-nearest-neighbors regression model for stability.
#' Differs from trainKNN (which does classification) by using continuous y.
#'
#' @param X Data frame of predictor features (rescaled env + coords)
#' @param y Numeric vector of stability scores
#' @return Trained caret model object
#' @keywords internal
#' @noRd
train_regression_knn <- function(X, y) {
  data_for_knn <- X
  data_for_knn$y <- y
  
  caret::train(
    y ~ .,
    data   = data_for_knn,
    method = "knn",
    tuneGrid = data.frame(k = c(3, 5, 7, 9)),
    trControl = caret::trainControl(method = "cv", number = 5)
  )
}


#' Extract top-3 most frequent cluster assignments per point across draws
#'
#' For each sample point, tallies which cluster it was assigned to across all
#' posterior draws, then extracts the 3 most frequent assignments and their
#' percentages.
#'
#' @param draw_clusterings integer matrix n_pts × n_draws
#' @param n_pts integer
#' @param n_draws integer
#' @return List with `labels_matrix` (n_pts × 3, cluster IDs) and
#'   `pct_matrix` (n_pts × 3, percentages)
#' @keywords internal
#' @noRd
compute_top3_clusters <- function(draw_clusterings, n_pts, n_draws) {
  labels_matrix <- matrix(NA_integer_, nrow = n_pts, ncol = 3)
  pct_matrix    <- matrix(NA_real_,    nrow = n_pts, ncol = 3)
  
  for (i in seq_len(n_pts)) {
    # Tally cluster assignments for point i across all draws
    assignments <- draw_clusterings[i, ]
    assignments <- assignments[!is.na(assignments)]
    
    if (length(assignments) == 0) next
    
    freq_table <- table(assignments)
    freq_table <- sort(freq_table, decreasing = TRUE)
    
    # Extract top 3 (or fewer if < 3 unique clusters were ever assigned)
    top_n <- min(3, length(freq_table))
    top_clusters <- as.integer(names(freq_table)[seq_len(top_n)])
    top_counts   <- as.numeric(freq_table[seq_len(top_n)])
    top_pcts     <- (top_counts / length(assignments)) * 100
    
    # Fill in the matrices
    labels_matrix[i, seq_len(top_n)] <- top_clusters
    pct_matrix[i, seq_len(top_n)]    <- top_pcts
    
    # If fewer than 3 unique clusters, replicate rank1 into empty slots
    if (top_n < 3) {
      for (j in (top_n + 1):3) {
        labels_matrix[i, j] <- top_clusters[1]
        pct_matrix[i, j]    <- top_pcts[1]
      }
    }
  }
  
  list(
    labels_matrix = labels_matrix,
    pct_matrix    = pct_matrix
  )
}



#'
#' @param model brmsfit
#' @param env_vars character vector of predictor names (no b_ prefix)
#' @param n_draws integer
#' @return matrix n_draws × length(env_vars), colnames = env_vars
#' @keywords internal
#' @noRd
extract_beta_draws <- function(model, env_vars, n_draws) {
  all_draws <- brms::as_draws_matrix(model)

  # brms stores fixed effects as "b_varname"
  b_cols <- paste0("b_", env_vars)
  found  <- b_cols[b_cols %in% colnames(all_draws)]

  if (length(found) == 0) {
    stop("No posterior draw columns matched environmental variable names.")
  }

  missing <- setdiff(b_cols, found)
  if (length(missing) > 0) {
    warning(sprintf(
      "These beta columns not found in draws: %s",
      paste(sub("^b_", "", missing), collapse = ", ")
    ))
  }

  n_total <- nrow(all_draws)
  idx     <- sample.int(n_total, size = min(n_draws, n_total), replace = FALSE)
  out     <- all_draws[idx, found, drop = FALSE]
  colnames(out) <- sub("^b_", "", colnames(out))
  as.matrix(out)
}


#' Rescale a point matrix by a specific draw's beta coefficients
#'
#' @param pt_env data frame of raw predictor values at sample points
#' @param env_vars character vector of variable names
#' @param var_mu named numeric vector of training means
#' @param var_sd named numeric vector of training SDs
#' @param betas named numeric vector of beta coefficients for this draw
#' @return data frame of scaled, beta-weighted values
#' @keywords internal
#' @noRd
rescale_points_by_betas <- function(pt_env, env_vars, var_mu, var_sd, betas) {
  out <- as.data.frame(matrix(
    NA_real_, nrow = nrow(pt_env), ncol = length(env_vars)
  ))
  colnames(out) <- env_vars

  for (v in env_vars) {
    out[[v]] <- ((pt_env[[v]] - var_mu[v]) / var_sd[v]) * abs(betas[v])
  }
  out
}


#' Add coordinate columns weighted to match the env layer range
#'
#' Mirrors the logic in [add_weighted_coordinates()] but operates on a
#' data frame of point values rather than a raster stack.
#'
#' @param sample_pts SpatVector of sample points
#' @param pt_env data frame of predictor values (no coordinate cols yet)
#' @param coord_wt numeric weight
#' @param env_vars character vector of variable names
#' @return data frame with additional columns coord_x_w and coord_y_w
#' @keywords internal
#' @noRd
add_coord_weights_to_points <- function(sample_pts, pt_env, coord_wt, env_vars) {
  coords  <- terra::crds(sample_pts)
  sdN     <- function(x) sqrt(mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))

  # Range of env variables (used as reference for weighting)
  env_ranges <- apply(pt_env[, env_vars, drop = FALSE], 2,
                      function(x) diff(range(x, na.rm = TRUE)))
  target_range <- max(env_ranges, na.rm = TRUE) * coord_wt

  # Scale each coordinate axis to target range
  scale_to_range <- function(v, tgt) {
    rng <- diff(range(v, na.rm = TRUE))
    if (rng == 0) return(rep(0, length(v)))
    (v - mean(v, na.rm = TRUE)) / rng * tgt
  }

  pt_env$coord_x_w <- scale_to_range(coords[, 1], target_range)
  pt_env$coord_y_w <- scale_to_range(coords[, 2], target_range)
  pt_env
}


#' Build symmetric co-occurrence probability matrix from draw clusterings
#'
#' For large n_pts this inner loop is the computational hotspot. The
#' implementation uses vectorised outer-equality checks per draw, which is
#' fast enough for n_pts ≤ 1000 and n_draws ≤ 200.
#'
#' @param draw_clusterings integer matrix n_pts × n_draws
#' @param n_pts integer
#' @param n_draws integer
#' @return numeric matrix n_pts × n_pts, values in [0, 1]
#' @keywords internal
#' @noRd
build_cooccurrence_matrix <- function(draw_clusterings, n_pts, n_draws) {
  co_mat <- matrix(0.0, nrow = n_pts, ncol = n_pts)

  for (d in seq_len(n_draws)) {
    labels   <- draw_clusterings[, d]
    # outer() produces an n_pts × n_pts logical matrix of same-cluster pairs
    same_clust <- outer(labels, labels, FUN = "==")
    same_clust[is.na(same_clust)] <- FALSE
    co_mat <- co_mat + same_clust
  }

  co_mat / n_draws
}


#' Compute per-point cluster stability scores
#'
#' Stability of point i = mean co-occurrence of i with all other points
#' sharing its consensus label. Points near cluster boundaries that were
#' frequently re-assigned in other draws will score low.
#'
#' @param co_mat numeric n_pts × n_pts co-occurrence matrix
#' @param consensus_labels integer vector length n_pts
#' @return numeric vector length n_pts, values in [0, 1]
#' @keywords internal
#' @noRd
compute_stability_scores <- function(co_mat, consensus_labels) {
  n_pts    <- nrow(co_mat)
  stability <- numeric(n_pts)

  for (i in seq_len(n_pts)) {
    cluster_i  <- consensus_labels[i]
    same_cluster <- which(consensus_labels == cluster_i)
    same_cluster <- setdiff(same_cluster, i)  # exclude self

    if (length(same_cluster) == 0) {
      stability[i] <- 1.0  # singleton: assign maximum stability
    } else {
      stability[i] <- mean(co_mat[i, same_cluster], na.rm = TRUE)
    }
  }
  stability
}


#' Interpolate point stability scores to a raster surface
#'
#' Uses inverse-distance-weighted interpolation (terra::interpIDW) to
#' produce a continuous stability surface from the discrete sample points.
#'
#' @param sample_pts SpatVector of sample points
#' @param stability_scores numeric vector, length = nrow(sample_pts)
#' @param mask_rast SpatRaster defining extent and resolution
#' @return SpatRaster named "cluster_stability", values in [0, 1]
#' @keywords internal
#' @noRd
interpolate_stability_to_raster <- function(sample_pts, stability_scores, mask_rast) {
  pts_df <- as.data.frame(terra::crds(sample_pts))
  pts_df$stability <- stability_scores

  stability_vect <- terra::vect(
    pts_df,
    geom  = c("x", "y"),
    crs   = terra::crs(mask_rast)
  )

  stab_rast <- terra::interpIDW(
    stability_vect,
    mask_rast,
    field = "stability",
    radius = terra::xres(mask_rast) * 10,
    power  = 2
  )
  stab_rast <- terra::mask(stab_rast, mask_rast)
  names(stab_rast) <- "cluster_stability"
  stab_rast
}


#' Project consensus clusters to full raster using KNN trained on original points
#'
#' Trains KNN classifiers using:
#' - Response: consensus cluster labels (rank1) and top-3 cluster rankings
#' - Features: posterior-mean-rescaled predictor values + weighted coordinates
#'
#' This ensures the spatial projection is faithful to the consensus clustering
#' rather than being a separate clustering exercise.
#'
#' @param sample_pts SpatVector of the original 500 fixed sample points
#' @param consensus_labels Integer vector of consensus cluster IDs (length = nrow(sample_pts))
#' @param stability_scores Numeric vector of stability scores (length = nrow(sample_pts))
#' @param top3_labels Integer matrix (n_pts × 3) of top-3 cluster assignments
#' @param predictors SpatRaster of raw environmental predictors
#' @param env_vars Character vector of predictor names
#' @param var_mu Named numeric vector of training means
#' @param var_sd Named numeric vector of training SDs
#' @param mean_betas Named numeric vector of posterior mean betas
#' @param coord_wt Numeric coordinate weight
#' @param mask_rast SpatRaster defining spatial domain
#' @param planar_proj CRS for planar projection
#' @return List with cluster_raster, stability_raster, rank2/3 rasters, knn models
#' @keywords internal
#' @noRd
project_consensus_to_raster <- function(
    sample_pts,
    consensus_labels,
    stability_scores,
    top3_labels,
    predictors,
    env_vars,
    var_mu,
    var_sd,
    mean_betas,
    coord_wt,
    mask_rast,
    planar_proj
) {
  # ── 1. Extract raw predictor values at the 500 fixed points ────────────────
  pt_env <- terra::extract(predictors[[env_vars]], sample_pts, ID = FALSE)
  
  # ── 2. Rescale using posterior mean betas ──────────────────────────────────
  pt_rescaled <- as.data.frame(matrix(
    NA_real_, nrow = nrow(pt_env), ncol = length(env_vars)
  ))
  colnames(pt_rescaled) <- env_vars
  
  for (v in env_vars) {
    pt_rescaled[[v]] <- ((pt_env[[v]] - var_mu[v]) / var_sd[v]) * abs(mean_betas[v])
  }
  
  # ── 3. Add weighted coordinates ────────────────────────────────────────────
  pt_rescaled <- add_coord_weights_to_points(
    sample_pts, pt_rescaled, coord_wt, env_vars
  )
  
  # ── 4. Attach response variables and filter complete cases ─────────────────
  pt_rescaled$rank1_cluster <- factor(consensus_labels)
  pt_rescaled$rank2_cluster <- factor(top3_labels[, 2])
  pt_rescaled$rank3_cluster <- factor(top3_labels[, 3])
  pt_rescaled$stability <- stability_scores
  
  complete_rows <- stats::complete.cases(pt_rescaled)
  pt_train <- pt_rescaled[complete_rows, ]
  
  if (nrow(pt_train) < 50) {
    stop("Fewer than 50 complete cases for KNN training. Check for NA values.")
  }
  
  # ── 5. Train KNNs ───────────────────────────────────────────────────────────
  cluster_features <- setdiff(names(pt_train), 
                               c("rank1_cluster", "rank2_cluster", "rank3_cluster", "stability"))
  
  # Rank1 (consensus) clusters
  knn_rank1 <- trainKNN(
    pt_train[, c(cluster_features, "rank1_cluster")],
    split_prop = 0.8
  )
  
  # Rank2 clusters (now always populated via replication)
  knn_rank2 <- trainKNN(
    pt_train[, c(cluster_features, "rank2_cluster")],
    split_prop = 0.8
  )
  
  # Rank3 clusters (now always populated via replication)
  knn_rank3 <- trainKNN(
    pt_train[, c(cluster_features, "rank3_cluster")],
    split_prop = 0.8
  )
  
  # Stability regression
  knn_stability <- train_regression_knn(
    X = pt_train[, cluster_features, drop = FALSE],
    y = pt_train$stability
  )
  
  # ── 6. Build posterior-mean-rescaled raster stack ──────────────────────────
  scaled_layers <- lapply(env_vars, function(v) {
    ((predictors[[v]] - var_mu[v]) / var_sd[v]) * abs(mean_betas[v])
  })
  pred_rescale <- terra::rast(scaled_layers)
  names(pred_rescale) <- env_vars
  
  # Add weighted coordinates to raster
  pred_rescale <- add_weighted_coordinates(pred_rescale, coord_wt)
  pred_rescale <- terra::mask(pred_rescale, mask_rast)
  
  # ── 7. Predict to full raster ──────────────────────────────────────────────
  cluster_raster <- terra::predict(
    pred_rescale, 
    model = knn_rank1$fit.knn, 
    na.rm = TRUE
  )
  
  rank2_raster <- terra::predict(
    pred_rescale, 
    model = knn_rank2$fit.knn, 
    na.rm = TRUE
  )
  names(rank2_raster) <- "rank2_cluster"
  
  rank3_raster <- terra::predict(
    pred_rescale, 
    model = knn_rank3$fit.knn, 
    na.rm = TRUE
  )
  names(rank3_raster) <- "rank3_cluster"
  
  stability_raster <- terra::predict(
    pred_rescale,
    model = knn_stability,
    na.rm = TRUE
  )
  names(stability_raster) <- "cluster_stability"
  
  # ── 8. Project to planar CRS ────────────────────────────────────────────────
  cluster_raster <- terra::project(cluster_raster, planar_proj)
  cluster_raster <- terra::mask(
    cluster_raster,
    terra::project(mask_rast, planar_proj)
  )
  
  stability_raster <- terra::project(stability_raster, planar_proj)
  stability_raster <- terra::mask(
    stability_raster,
    terra::project(mask_rast, planar_proj)
  )
  
  rank2_raster <- terra::project(rank2_raster, planar_proj)
  rank2_raster <- terra::mask(rank2_raster, terra::project(mask_rast, planar_proj))
  
  rank3_raster <- terra::project(rank3_raster, planar_proj)
  rank3_raster <- terra::mask(rank3_raster, terra::project(mask_rast, planar_proj))
  
  list(
    cluster_raster   = cluster_raster,
    stability_raster = stability_raster,
    rank2_raster     = rank2_raster,
    rank3_raster     = rank3_raster,
    knn_cluster      = knn_rank1$fit.knn,
    knn_rank2        = knn_rank2$fit.knn,
    knn_rank3        = knn_rank3$fit.knn,
    knn_stability    = knn_stability
  )
}