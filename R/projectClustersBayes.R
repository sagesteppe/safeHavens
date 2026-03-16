#' Project posterior-clustered environmental space onto a future climate scenario
#'
#' Bayesian analogue of [projectClusters()]. The steps are:
#'
#' 1. Rescale future predictors via `RescaleRasters_bayes` to match the feature space
#'    that [PosteriorCluster()] used when training the consensus KNN.
#' 2. Run `dismo::mess` to identify novel climate cells.
#' 3. Predict known-climate cells with the consensus KNN classifier from
#'    [PosteriorCluster()].
#' 4. If novel cells exist and `cluster_novel = TRUE`, cluster them
#'    independently with `NbClust` / KNN (`cluster_novel_areas`).
#' 5. Sample points from every cluster (existing + novel), build a single
#'    tree, and extract nearest-existing-cluster relationships via silhouette
#'    (`analyze_cluster_relationships`).
#' 6. Project future stability by replaying posterior beta draws on the future
#'    climate surface (`project_future_draws`). Stability here is a genuinely
#'    forward-looking quantity: the fraction of draws in which each cell's
#'    modal cluster assignment was agreed upon. This sidesteps the concern of
#'    transferring a current-era stability surface trained in current climate
#'    space onto a future surface that may differ substantially.
#' 7. Polygonise and calculate area / centroid changes.
#'
#' @param bSDM_object Output list from [bayesianSDM()].
#' @param posterior_clusters Output list from [PosteriorCluster()].
#' @param future_predictors `SpatRaster` of future climate.
#'   Layer names must match those used in `bSDM_object`.
#' @param current_predictors `SpatRaster` of current climate (used for MESS
#'   reference and MESS computation (raw units required).
#' @param threshold_rasts Output list from [PostProcessSDM()], containing:
#'   `$FinalRasters` (a `SpatRaster` stack including a binary threshold layer)
#'   and `$Threshold` (a one-row `data.frame` of threshold metric values).
#' @param planar_proj EPSG code or proj4 string for a planar projection in
#'   metres (same as used in the current analysis).
#' @param coord_wt Ignored. Coordinate weighting is read directly from
#'   `posterior_clusters$ScalingParams$coord_wt` to guarantee the future
#'   feature space matches the scale the KNNs were trained on. The parameter
#'   is retained in the signature for backwards compatibility only.
#' @param mess_threshold MESS values below this are treated as novel climate.
#'   Default `0`.
#' @param cluster_novel Logical, default `TRUE`. If `TRUE` and novel cells
#'   exist, cluster them independently with NbClust + KNN. If `FALSE` novel
#'   cells are left as `NA`.
#' @param n_novel_pts Numeric, default `500`. Number of points to sample from
#'   novel areas for independent clustering.
#' @param n_sample_per_cluster Number of points to sample from each cluster
#'   (existing *and* novel) for the silhouette relationship analysis.
#'   Default `50`.
#' @param n_future_draws Integer. Number of posterior beta draws to use when
#'   computing future cluster-assignment stability. Defaults to the number of
#'   draws stored in `posterior_clusters$ScalingParams$beta_draws` (i.e.
#'   however many were used in [PosteriorCluster()]), capped at `100`. Pass
#'   an explicit integer to override.
#' @param n_future_pts Integer. Number of points to sample from the future
#'   suitable-habitat surface for the draw loop. These are analogous to the
#'   fixed sample points in [PosteriorCluster()] and serve the same purpose:
#'   stability is computed on this point set, then painted to raster once via
#'   KNN. Default `500`.
#' @param nbclust_args Named list of arguments forwarded to
#'   `NbClust::NbClust`. Sensible defaults are set internally (`min.nc = 2`,
#'   `max.nc = 10`, `method = "ward.D2"`, `index = "all"`).
#' @param thresh_metric Character. Default `"sensitivity"`. The
#'   `dismo::threshold` metric used to cut the future SDM into a binary
#'   surface.
#'
#' @returns A named list:
#'   \item{`clusters_raster`}{Smoothed `SpatRaster` of consensus cluster IDs.}
#'   \item{`clusters_sf`}{`sf` polygons with column `ID`.}
#'   \item{`suitable_habitat`}{`SpatRaster` - binary suitable habitat under
#'     future conditions.}
#'   \item{`novel_mask`}{`SpatRaster` - `1` where MESS < threshold.}
#'   \item{`mess`}{`SpatRaster` - raw MESS scores.}
#'   \item{`stability`}{`SpatRaster` - draw-based cluster-assignment stability
#'     under future climate. Fraction of posterior draws agreeing on each
#'     cell's modal cluster. Values in (0, 1); low values flag cells where
#'     beta uncertainty produces ambiguous cluster assignment in future space.
#'     Unlike the current-era stability surface, this is computed directly in
#'     future climate space and is not a transfer from [PosteriorCluster()].}
#'   \item{`changes`}{`data.frame` - per-cluster area and centroid-shift
#'     metrics.}
#'   \item{`novel_similarity`}{`data.frame` - nearest existing cluster and
#'     mean silhouette width for each novel cluster. Zero-row if none.}
#'
#' @seealso [projectClusters()] for the elastic-net equivalent,
#'   [PosteriorCluster()], [bayesianSDM()]
#' @export
projectClustersBayes <- function(
  bSDM_object,
  posterior_clusters,
  future_predictors,
  current_predictors,
  threshold_rasts,
  planar_proj,
  coord_wt             = NULL,   # ignored; read from ScalingParams
  mess_threshold       = 0,
  cluster_novel        = TRUE,
  n_novel_pts          = 500,
  n_sample_per_cluster = 50,
  n_future_draws       = NULL,
  n_future_pts         = 500,
  nbclust_args         = list(),
  thresh_metric        = "sensitivity"
) {
  # --- Unpack bSDM_object ------------------------------------------------------
  current_model <- bSDM_object$Model
  pred_mat      <- bSDM_object$PredictMatrix

  # --- Unpack posterior_clusters ------------------------------------------------
  knn_consensus <- posterior_clusters$KNNModels$KNN_Cluster
  scaling       <- posterior_clusters$ScalingParams
  mean_betas    <- scaling$mean_betas
  beta_draws    <- scaling$beta_draws   # n_stored_draws x n_vars matrix
  coord_wt      <- scaling$coord_wt     # must match what KNNs were trained on
  env_vars      <- names(mean_betas)

  # Resolve n_future_draws: default to however many draws came in, cap at 100
  n_stored       <- nrow(beta_draws)
  n_future_draws <- if (is.null(n_future_draws)) {
    min(n_stored, 100L)
  } else {
    min(as.integer(n_future_draws), n_stored)
  }

  # --- 1. SDM prediction on future climate (population-level, no GP) ------------
  # The spatial smooth (gp_x, gp_y) is kept active (re_formula = NULL) by
  # supplying proper planar-km coordinate rasters built from the future extent.
  # This mirrors create_bayes_spatial_predictions exactly. The GP extrapolates
  # into future space, which is imperfect but far more honest than zeroing
  # coordinates. The threshold was calibrated on the GP-inclusive current model
  # so this keeps the prediction on a consistent scale.
  thresholds <- threshold_rasts$Threshold

  # Build gp_x / gp_y coordinate rasters for the future extent in planar km,
  # matching exactly what create_bayes_spatial_predictions does for current data.
  # re_formula = NULL (default) keeps the spatial smooth active - extrapolating
  # the GP into future space is more honest than zeroing coordinates.
  future_pred_stack <- .build_future_pred_stack(
    future_predictors = future_predictors,
    env_vars          = env_vars,
    planar_proj       = planar_proj
  )

  n_iter <- nrow(beta_draws)   #
  predict_fun <- .make_brms_predict_fun(n_iter)

  future_sdm_stack <- terra::predict(
    future_pred_stack,
    model = current_model,
    fun   = predict_fun,
    na.rm = TRUE,
    cores = 1
  )
  future_sdm <- future_sdm_stack[[1]]  # mean layer

  cut              <- thresholds[[thresh_metric]]
  suitable_habitat <- future_sdm >= cut

  # --- 2. MESS - identify novel climate --------------------------------------------
  ref_matrix <- terra::as.data.frame(
    current_predictors[[env_vars]],
    na.rm = TRUE
  )

  mess_result <- terra::rast(
    dismo::mess(
      x    = raster::stack(future_predictors[[env_vars]]),
      v    = ref_matrix,
      full = FALSE
    )
  )
  mess_overall  <- mess_result[[terra::nlyr(mess_result)]]
  novel_climate <- mess_overall < mess_threshold

  # --- 3. Partition suitable habitat ------------------------------------------
  suitable_known          <- terra::ifel(suitable_habitat & !novel_climate, 1, NA)
  suitable_novel          <- terra::ifel(suitable_habitat &  novel_climate, 1, NA)
  suitable_habitat_binary <- terra::ifel(suitable_habitat, 1, NA)
  novel_mask_binary       <- terra::ifel(novel_climate,    1, NA)

  # --- 4. Rescale future predictors to match current-era KNN feature space -------
  # PosteriorCluster received predictors already processed by RescaleRasters_bayes.
  # The KNNs were therefore trained on that rescaled space. We apply the same
  # transformation to future_predictors here so the feature space matches.
  # rescaleFutureBayes (which used raw pred_mat moments) is intentionally not
  # used - it produced the wrong scale because it was unaware of the upstream
  # RescaleRasters_bayes transformation.
  # nocov start
  future_rr <- RescaleRasters_bayes(
    model         = current_model,
    predictors    = future_predictors,
    training_data = bSDM_object$TrainData,
    pred_mat      = pred_mat
  )
  future_rescaled <- future_rr$RescaledPredictors
  future_rescaled <- add_weighted_coordinates(future_rescaled, coord_wt)
  names(future_rescaled)[names(future_rescaled) == "x"] <- "coord_x_w"
  names(future_rescaled)[names(future_rescaled) == "y"] <- "coord_y_w"
  # nocov end

  # --- 5. Project consensus clusters onto suitable + known areas -----------
  known_clusters <- terra::predict(
    future_rescaled,
    model = knn_consensus,
    na.rm = TRUE
  )
  known_clusters <- terra::mask(known_clusters, suitable_known)

  n_novel_cells <- terra::global(suitable_novel, "sum", na.rm = TRUE)[[1]]

  # --- 6. Cluster novel areas (NbClust + KNN, same as enet path) -----------
  if (cluster_novel && !is.na(n_novel_cells) && n_novel_cells > 0) {
    # nocov start
    novel_result <- cluster_novel_areas(
      future_rescaled = future_rescaled,
      novel_mask      = suitable_novel,
      n_novel_pts     = n_novel_pts,
      next_cluster_id = max(posterior_clusters$Geometry$ID) + 1,
      nbclust_args    = nbclust_args
    )
    # nocov end

    # novel_clusters already carries offset IDs (next_cluster_id was set to
    # max(Geometry$ID) + 1 in the cluster_novel_areas call above) - do NOT
    # add the offset again here, matching the enet path in projectClusters.
    novel_clusters <- novel_result$clusters_raster
    known_numeric  <- terra::as.int(known_clusters)
    novel_numeric  <- terra::as.int(novel_clusters)
    final_clusters <- terra::cover(known_numeric, novel_numeric)

    # --- 7. Relationship analysis ----------------------------------------------
    # nocov start
    novel_similarity <- analyze_cluster_relationships(
      clusters_raster      = final_clusters,
      future_rescaled      = future_rescaled,
      existing_ids         = unique(posterior_clusters$Geometry$ID),
      novel_ids            = as.numeric(unlist(terra::unique(novel_numeric))),
      n_sample_per_cluster = n_sample_per_cluster
    )
    # nocov end

  } else {
    final_clusters   <- known_clusters
    novel_similarity <- data.frame(
      novel_cluster_id     = integer(),
      nearest_existing_id  = integer(),
      avg_silhouette_width = numeric()
    )
  }

  # -- 8. Draw-based future stability ---------------------------------------------
  # Replay posterior beta draws on the future climate surface. Stability at
  # each cell = fraction of draws agreeing on the modal cluster assignment.
  # Computed entirely in future climate space - cells in novel or strongly
  # shifted climate will naturally show low draw-agreement rather than
  # inheriting a stability score from the current era.
  # nocov start
  message(sprintf(
    "Computing future stability over %d posterior draws ...", n_future_draws
  ))
  future_stability <- project_future_draws(
    future_rescaled_rr   = future_rr,
    suitable_mask        = suitable_habitat_binary,
    env_vars             = env_vars,
    beta_draws           = beta_draws[seq_len(n_future_draws), , drop = FALSE],
    coord_wt             = coord_wt,
    knn_consensus        = knn_consensus,
    n_future_pts         = n_future_pts,
    existing_cluster_ids = unique(posterior_clusters$Geometry$ID)
  )
  # nocov end

  # --- 9. Smooth noise ----------------------------------------------------------
  final_clean <- final_clusters

  # --- 10. Polygonize ----------------------------------------------------------
  # IDs are kept as-is: current-era clusters retain their geographic IDs from
  # PosteriorCluster; novel clusters carry offset IDs (max(current) + 1 ...).
  # calculate_changes() matches by ID so lost/gained clusters are detected
  # correctly without any renumbering - identical to the enet path.
  clusters_sf <- terra::as.polygons(final_clean) |>
    sf::st_as_sf() |>
    sf::st_make_valid()
  names(clusters_sf)[1] <- "ID"
  clusters_sf$ID <- as.integer(clusters_sf$ID)

  # --- 11. Change metrics ----------------------------------------------------------
  # nocov start
  changes <- calculate_changes(
    current_sf  = posterior_clusters$Geometry,
    future_sf   = clusters_sf,
    planar_proj = planar_proj
  )
  # nocov end

  # --- 12. Return ----------------------------------------------------------
  list(
    clusters_raster  = final_clean,
    Geometry         = clusters_sf,
    suitable_habitat = suitable_habitat_binary,
    novel_mask       = novel_mask_binary,
    mess             = mess_overall,
    stability        = future_stability,
    changes          = changes,
    novel_similarity = novel_similarity
  )
}


# ==============================================================================
#  Draw-based future stability
# ==============================================================================

#' Compute draw-based cluster-assignment stability in future climate space
#'
#' Samples a fixed point set from the future suitable-habitat surface, then
#' for each posterior beta draw reweights the already-RescaleRasters_bayes-scaled
#' point values by the draw's betas (as a ratio over the mean beta) and predicts
#' cluster assignment via the consensus KNN. Stability at each point is the
#' fraction of draws that agreed on that point's modal cluster.
#'
#' Critically, point values are extracted from `future_rescaled_rr$RescaledPredictors`
#' - the output of RescaleRasters_bayes applied to future climate - so the feature
#' space is on exactly the same scale the KNN was trained on. The per-draw
#' reweighting then modulates those values by the ratio of the draw's beta to the
#' posterior mean beta, mirroring what PosteriorCluster does internally.
#'
#' @param future_rescaled_rr Return value of `RescaleRasters_bayes()` applied to
#'   future predictors. Must contain `$RescaledPredictors` (a `SpatRaster`).
#' @param suitable_mask `SpatRaster` - binary suitable habitat (1 / NA).
#' @param env_vars Character vector of environmental predictor names.
#' @param beta_draws Numeric matrix (n_draws x n_vars) of posterior beta samples.
#' @param coord_wt Numeric coordinate weight matching [PosteriorCluster()].
#' @param knn_consensus Trained caret KNN (`KNN_Cluster` from [PosteriorCluster()]).
#' @param n_future_pts Integer. Points to sample from future surface. Default `500`.
#' @param existing_cluster_ids Integer vector of valid current-era cluster IDs.
#'
#' @returns `SpatRaster` named `"future_stability"`, values in (0, 1).
#' @keywords internal
#' @noRd
project_future_draws <- function(
  future_rescaled_rr,
  suitable_mask,
  env_vars,
  beta_draws,
  coord_wt,
  knn_consensus,
  n_future_pts         = 500L,
  existing_cluster_ids
) {
  n_draws         <- nrow(beta_draws)
  future_rescaled <- future_rescaled_rr$RescaledPredictors
  mean_betas_vec  <- colMeans(beta_draws)

  # --- 1. Sample fixed points from future suitable habitat ------------------
  message("  Sampling fixed future point set ...")

  n_available <- terra::global(suitable_mask, "sum", na.rm = TRUE)[[1]]
  if (is.na(n_available) || n_available < 10) {
    warning(
      "project_future_draws: no suitable-habitat cells found. ",
      "Returning NA stability raster."
    )
    out <- suitable_mask * NA_real_
    names(out) <- "future_stability"
    return(out)
  }

  n_sample <- min(n_future_pts, floor(n_available))
  if (n_sample < n_future_pts) {
    message(sprintf(
      "  Note: only %d suitable cells available; clamping n_future_pts from %d to %d.",
      n_sample, n_future_pts, n_sample
    ))
  }

  sample_pts <- terra::spatSample(
    suitable_mask,
    size      = n_sample,
    method    = "random",
    as.points = TRUE,
    na.rm     = TRUE
  )

  # Extract already-RescaleRasters_bayes-scaled values - same space as KNN training
  pt_rr         <- terra::extract(future_rescaled[[env_vars]], sample_pts, ID = FALSE)
  complete_rows <- stats::complete.cases(pt_rr)
  sample_pts    <- sample_pts[complete_rows, ]
  pt_rr         <- pt_rr[complete_rows, , drop = FALSE]
  n_pts_actual  <- nrow(pt_rr)

  if (n_pts_actual < 10L) {
    warning(
      "project_future_draws: fewer than 10 points with complete data. ",
      "Returning NA stability raster."
    )
    out <- suitable_mask * NA_real_
    names(out) <- "future_stability"
    return(out)
  }

  message(sprintf(
    "  %d of %d available future points have complete predictor data.",
    n_pts_actual, n_sample
  ))

  # --- 2. Draw loop ------------------------------------------------------
  # Per-draw reweighting: the raster values are already scaled by |mean_beta|
  # via RescaleRasters_bayes. For draw d we rescale by |beta_d| instead, which
  # is equivalent to multiplying by |beta_d| / |mean_beta| per variable.
  # This mirrors exactly what PosteriorCluster does in its draw loop via
  # rescale_points_by_betas, but starting from the correct feature space.
  draw_assignments <- matrix(NA_integer_, nrow = n_pts_actual, ncol = n_draws)

  for (d in seq_len(n_draws)) {
    if (d %% 25 == 0) message(sprintf("  Draw %d / %d", d, n_draws))

    betas_d <- stats::setNames(
      as.numeric(beta_draws[d, , drop = FALSE]),
      colnames(beta_draws)
    )

    # Reweight: divide out mean beta, multiply in draw beta
    pt_scaled <- pt_rr
    for (v in env_vars) {
      mb <- mean_betas_vec[v]
      pt_scaled[[v]] <- if (abs(mb) < .Machine$double.eps) {
        pt_rr[[v]] * abs(betas_d[v])
      } else {
        pt_rr[[v]] * (abs(betas_d[v]) / abs(mb))
      }
    }

    pt_full <- add_coord_weights_to_points(
      sample_pts = sample_pts,
      pt_env     = pt_scaled,
      coord_wt   = coord_wt,
      env_vars   = env_vars
    )
    pt_full <- pt_full[stats::complete.cases(pt_full), , drop = FALSE]

    preds     <- stats::predict(knn_consensus, newdata = pt_full)
    preds_int <- suppressWarnings(as.integer(as.character(preds)))
    preds_int[!preds_int %in% existing_cluster_ids] <- NA_integer_
    draw_assignments[seq_len(nrow(pt_full)), d] <- preds_int
  }

  # --- 3. Modal cluster and per-point stability ---------------------------------
  modal_cluster <- apply(draw_assignments, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) return(NA_integer_)
    as.integer(names(sort(table(x), decreasing = TRUE)[1L]))
  })

  stability_scores <- vapply(seq_len(n_pts_actual), function(i) {
    x     <- draw_assignments[i, ]
    x     <- x[!is.na(x)]
    modal <- modal_cluster[i]
    if (length(x) == 0L || is.na(modal)) return(NA_real_)
    mean(x == modal)
  }, numeric(1L))

  # --- 4. Paint stability to raster - one terra::predict call ------------------
  # Use the already-correctly-scaled future_rescaled + coordinates as features,
  # mirroring exactly the feature space the stability regression will generalise over.
  pt_mean_scaled <- pt_rr
  pt_mean_scaled <- add_coord_weights_to_points(
    sample_pts = sample_pts,
    pt_env     = pt_mean_scaled,
    coord_wt   = coord_wt,
    env_vars   = env_vars
  )

  complete_train <- stats::complete.cases(pt_mean_scaled) & !is.na(stability_scores)
  pt_train       <- pt_mean_scaled[complete_train, , drop = FALSE]
  stab_train     <- stability_scores[complete_train]

  if (nrow(pt_train) < 10L) {
    warning(
      "project_future_draws: too few complete training points for stability KNN. ",
      "Returning NA stability raster."
    )
    out <- suitable_mask * NA_real_
    names(out) <- "future_stability"
    return(out)
  }

  knn_future_stab <- train_regression_knn(X = pt_train, y = stab_train)

  # Prediction raster: future_rescaled (already correct scale) + coord layers
  pred_rescale_future <- future_rescaled[[env_vars]]
  pred_rescale_future <- add_weighted_coordinates(pred_rescale_future, coord_wt)
  names(pred_rescale_future)[names(pred_rescale_future) == "x"] <- "coord_x_w"
  names(pred_rescale_future)[names(pred_rescale_future) == "y"] <- "coord_y_w"
  pred_rescale_future <- terra::mask(pred_rescale_future, suitable_mask)

  stability_raster <- terra::predict(
    pred_rescale_future,
    model = knn_future_stab,
    na.rm = TRUE
  )
  stability_raster <- terra::clamp(stability_raster, lower = 0, upper = 1)
  names(stability_raster) <- "future_stability"
  stability_raster
}

# ==============================================================================
#  brms prediction helper
# ==============================================================================

#' Build a future predictor stack with gp_x / gp_y coordinate layers
#'
#' Mirrors the coordinate-raster construction in `create_bayes_spatial_predictions`
#' so that the brms spatial smooth receives properly projected km-scale
#' coordinates rather than dummy zeros. Must be called with the same
#' `planar_proj` used when fitting the model.
#'
#' @param future_predictors `SpatRaster` of future climate.
#' @param env_vars Character vector of env predictor names to include.
#' @param planar_proj EPSG or proj4 string for planar km projection.
#' @returns `SpatRaster` with env layers + `gp_x` + `gp_y`.
#' @keywords internal
#' @noRd
.build_future_pred_stack <- function(future_predictors, env_vars, planar_proj) {
  template    <- future_predictors[[env_vars[[1]]]]
  coords_ll   <- terra::as.data.frame(template, xy = TRUE, cells = FALSE)[, c("x", "y")]

  coords_vect <- terra::vect(
    coords_ll,
    geom = c("x", "y"),
    crs  = terra::crs(future_predictors)
  )
  coords_proj <- terra::project(coords_vect, planar_proj)
  coords_km   <- terra::crds(coords_proj) / 1000

  gp_x_rast <- terra::rast(template)
  gp_y_rast <- terra::rast(template)
  terra::values(gp_x_rast) <- coords_km[, 1]
  terra::values(gp_y_rast) <- coords_km[, 2]
  names(gp_x_rast) <- "gp_x"
  names(gp_y_rast) <- "gp_y"

  c(future_predictors[[env_vars]], gp_x_rast, gp_y_rast)
}


#' Return a brms mean-prediction function closed over n_iter
#'
#' Returns a function compatible with `terra::predict(..., fun = ...)` that
#' computes posterior mean (and SD) predictions. Using a closure lets us
#' capture `n_iter` without needing an extra argument channel that
#' `terra::predict` does not support.
#'
#' @param n_iter Total posterior draws available; half are used for speed,
#'   matching the convention in `create_bayes_spatial_predictions`.
#' @returns A function(model, data, ...) returning a 2-column matrix
#'   (mean, sd).
#' @keywords internal
#' @noRd
.make_brms_predict_fun <- function(n_iter) {
  function(model, data, ...) {
    epred <- brms::posterior_epred(
      model,
      newdata          = data,
      allow_new_levels = TRUE,
      ndraws           = round(n_iter / 2, 0)
    )
    cbind(
      mean = colMeans(epred),
      sd   = apply(epred, 2, stats::sd)
    )
  }
}