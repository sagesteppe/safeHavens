#' Project current environmental clusters onto a future climate scenario
#'
#' Main entry point for the future-projection workflow.  The steps are:
#'
#' 1. Rescale future predictors with current betas (`rescaleFuture`).
#' 2. Run `dismo::mess` to identify novel climate cells.
#' 3. Predict known-climate cells with the existing KNN classifier.
#' 4. If novel cells exist and `cluster_novel = TRUE`, cluster them
#'    independently with `NbClust` (`cluster_novel_areas`).
#' 5. Sample points from every cluster (existing + novel), force them onto a
#'    single tree, and extract nearest-existing-cluster relationships via
#'    silhouette (`analyze_cluster_relationships`).
#' 6. Polygonise and calculate area / centroid changes.
#'
#' @param eSDM_object output from `elasticSDM()`.
#' @param current_clusters Output list from `EnvironmentalBasedSample()`. 
#' @param future_predictors SpatRaster of future climate.  
#' Layer names must match those retained in `current_model`.
#' @param current_predictors SpatRaster of current climate (used for MESS reference and for `rescaleFuture` standardisation).
#' @param planar_proj EPSG code or proj4 string for a planar projection in metres (same as used in the current analysis).
#' @param coord_wt Numeric, default `0.001`. Coordinate weighting passed to `add_weighted_coordinates`. 
#' @param mess_threshold MESS values below this are treated as novel climate. Default `0`.
#' @param cluster_novel Boolean, defualt `TRUE`. If `TRUE` and novel cells exist, cluster them independently with NbClust.  
#'  If `FALSE` novel cells are left as `NA`.
#' @param n_novel_pts Numeric, default 500. Number of points to sample from novel areas for clustering. 
#' @param n_sample_per_cluster Number of points to sample from each cluster
#'   (existing *and* novel) for the relationship tree.  Default `50`.
#' @param nbclust_args Named list of arguments forwarded to `NbClust::NbClust`.  
#' Sensible defaults are set internally (`min.nc = 2`, `max.nc = 10`, `method = "ward.D2"`, `index = "all"`).
#'
#' @returns A named list:
#'   \item{clusters_sf}{`sf` polygons with column `ID`.}
#'   \item{suitable_habitat}{Raster of masked suitable habitat under future conditions}
#'   \item{novel_mask}{`SpatRaster` — logical, `TRUE` where MESS < threshold.}
#'   \item{mess}{`SpatRaster` — raw MESS scores (minimum across all variables).}
#'   \item{changes}{`tibble` — per-cluster area and centroid-shift metrics.}
#'   \item{novel_similarity}{`tibble` — nearest existing cluster and mean
#'     silhouette width for each novel cluster.  Zero-row tibble if none.}
#'
#' @export
projectClusters <- function(
  eSDM_object,
  current_clusters,
  future_predictors,
  current_predictors,
  planar_proj,
  coord_wt            = 0.001,
  mess_threshold      = 0,
  cluster_novel       = TRUE,
  n_novel_pts         = 500,
  n_sample_per_cluster = 50,
  nbclust_args        = list(),
  thresh_metric       = 'sensitivity',
  thresholds
) {

  ## subset some items from eSDMmodel
  test_points = eSDM_object$TestData
  train_points = eSDM_object$TrainData
  predict_matrix = eSDM_object$PredictMatrix
  current_model = eSDM_object$Model

  # ---- 0. Extract model variables ----
  coef_mat <- as.matrix(stats::coef(current_model))
  vars <- rownames(coef_mat)#[coef_mat[, 1] != 0]
  vars <- vars[vars != "(Intercept)"]

  # ---- 1. SDM prediction on ORIGINAL future climate ----
  predfun <- function(model, data, ...) {
    stats::predict(model, newx = as.matrix(data), type = 'response')
  }

  future_sdm <- terra::predict(
    future_predictors[[vars]], 
    model = current_model, 
    fun = predfun, 
    na.rm = TRUE
  )

  # Threshold to binary habitat suitability
  cut <- thresholds[[thresh_metric]] 
  suitable_habitat <- future_sdm >= cut

  # ---- 2. MESS — identify novel climate (species-relevant vars only) ----
  ref_matrix <- terra::as.data.frame(
    current_predictors[[vars]], 
    na.rm = TRUE
  )

  mess_result <- terra::rast(  
    dismo::mess(
      x = raster::stack(future_predictors[[vars]]),
      v = ref_matrix,
      full = FALSE
    )
  )
  mess_overall <- mess_result[[terra::nlyr(mess_result)]] 
  novel_climate <- mess_overall < mess_threshold

  # ---- 3. similar habitat to existing, but forecast ----
  suitable_known <- terra::ifel(suitable_habitat & !novel_climate, 1, NA)

  # Suitable habitat that is in NOVEL climate  
  suitable_novel <- terra::ifel(suitable_habitat & novel_climate, 1, NA)

  # Binary versions for output
  suitable_habitat_binary <- terra::ifel(suitable_habitat, 1, NA)
  novel_mask_binary <- terra::ifel(novel_climate, 1, NA)
  
  # ---- 4. Rescale future for clustering (KNN needs rescaled data) ----
  future_rescaled <- rescaleFuture(
    model = current_model,
    future_predictors = future_predictors,
    current_predictors = current_predictors, 
    training_data = train_points,
    pred_mat = predict_matrix
  )
  
  # Add weighted coordinates
  future_rescaled <- add_weighted_coordinates(future_rescaled, coord_wt)

  # ---- 5. Project t0 clusters into suitable+known areas ----
  known_clusters <- terra::predict(
    future_rescaled,
    model = current_clusters$fit.knn,
    na.rm = TRUE
  )  

  # Restrict to suitable+known zone only
  known_clusters <- terra::mask(known_clusters, suitable_known)
  
  # Restrict to suitable+known zone only
  areas_needing_clustering <- terra::ifel(
    !is.na(suitable_novel) & is.na(known_clusters),
    1,
    NA
  )
  n_novel_cells <- terra::global(suitable_novel, "sum", na.rm = TRUE)[[1]]
  
  # ---- 6. Cluster suitable+novel areas independently ----
  suitable_habitat_binary <- terra::ifel(suitable_habitat, 1, NA)

  if (cluster_novel && n_novel_cells > 0) {

    novel_result <- cluster_novel_areas(
      future_rescaled  = future_rescaled,
      novel_mask       = suitable_novel,  # Only suitable novel areas
      n_novel_pts      = n_novel_pts,
      next_cluster_id  = max(current_clusters$Geometry$ID) + 1,
      nbclust_args     = nbclust_args
    )

    novel_clusters <- novel_result$clusters_raster
    # ---- Merge: combine known + novel ----
    
    # Overlay novel clusters (they should be in different areas due to masks)
    # Use terra::cover to fill NAs in known_clusters with novel values
    known_numeric <- terra::as.int(known_clusters)
    novel_numeric <- terra::as.int(novel_clusters)
    novel_numeric <- novel_numeric + max(current_clusters$Geometry$ID)
    final_clusters <- terra::cover(known_numeric, novel_numeric)

    # ---- 7. Relationship analysis ----
    novel_similarity <- analyze_cluster_relationships(
      clusters_raster      = final_clusters,
      future_rescaled      = future_rescaled,
      existing_ids         = unique(current_clusters$Geometry$ID),
      novel_ids            = as.numeric(unlist(unique(novel_numeric))),
      n_sample_per_cluster = n_sample_per_cluster
    )

  } else {
    final_clusters   <- known_clusters
    novel_similarity <- tibble::tibble(
      novel_cluster_id     = integer(),
      nearest_existing_id  = integer(),
      avg_silhouette_width = numeric()
    )
  }

  # ---- 8. Smooth noise ----
  template <- terra::rast(final_clusters)
  final_clusters <- mask(final_clusters, suitable_habitat_binary)
  agg <- terra::aggregate(
    final_clusters,
    fact = 2,
    fun  = modal,
    na.rm = TRUE
  )

  final_clean <- terra::resample(agg, template, method = "near")
  
  # ---- 9. Polygonize ----
  clusters_sf <- terra::as.polygons(final_clean) |>
    sf::st_as_sf() |>
    sf::st_make_valid()
  names(clusters_sf)[1] <- "ID"

  # ---- 10. Change metrics ----
  changes <- calculate_changes(
    current_sf  = current_clusters$Geometry,
    future_sf   = clusters_sf,
    planar_proj = planar_proj
  )

  # ---- 11. Return ----
  list(
    clusters_raster  = final_clean,
    clusters_sf      = clusters_sf,
    suitable_habitat = suitable_habitat_binary,
    novel_mask       = novel_mask_binary,
    mess             = mess_overall,
    changes          = changes,
    novel_similarity = novel_similarity
  )
}
# ------------------------------------------------------------------------------

#' Rescale future climate predictors using current model coefficients
#'
#' Standardises future predictors using the mean and SD from *current* climate,
#' then applies the glmnet beta weights.  This is the same weighting applied
#' upstream in `RescaleRasters` — here we just do it against a different
#' raster stack.
#'
#' @param model glmnet model object from `elasticSDM_noPCNM()$Model`.
#' @param future_predictors SpatRaster of future climate variables.  Names must
#'   match those retained in `model`.
#' @param current_predictors SpatRaster of current climate (provides the
#'   standardisation parameters).
#' @param training_data the same data that went into the glmnet model, this is used
#' for calculating variance which is required for the scaling process. From `elasticSDM`
#' @param pred_mat the Prediction matrix from `elasticSDM`
#' 
#' @returns SpatRaster with rescaled future predictors (one layer per
#'   non-zero coefficient).
#' @export
rescaleFuture <- function(model, future_predictors, current_predictors, training_data, pred_mat) {
  
  # Custom SD function (matches RescaleRasters)
  sdN <- function(x) {
    sqrt((1 / length(x)) * sum((x - mean(x))^2))
  }
  
  # ---- 1. Extract coefficients (same as RescaleRasters) ----
  coef_mat <- as.matrix(stats::coef(model))
  coef_vec <- coef_mat[, 1]
  coef_names <- rownames(coef_mat)
  coef_tab <- data.frame(
    Variable = coef_names,
    Coefficient = coef_vec,
    row.names = NULL
  )
  
  # Drop intercept
  coef_tab <- coef_tab[coef_tab$Variable != "(Intercept)", , drop = FALSE]
  
  # ---- 2. Response variance (from current/training) ----
  yvar <- sdN(as.numeric(training_data$occurrence) - 1)
  
  # ---- 3. Predictor SDs from CURRENT pred_mat ----
  stopifnot(all(coef_tab$Variable %in% colnames(pred_mat)))
  
  x_sd <- vapply(
    coef_tab$Variable,
    function(v) sdN(pred_mat[, v]),
    numeric(1)
  )
  
  # ---- 4. Compute beta coefficients ----
  coef_tab$BetaCoefficient <- (coef_tab$Coefficient / yvar) * x_sd
  
  # ---- 5. Rescale FUTURE predictors using CURRENT parameters ----
  future_sub <- terra::subset(future_predictors, coef_tab$Variable)
  layer_names <- names(future_sub)
  
  out_layers <- lapply(layer_names, function(lyr_name) {
    # Use CURRENT pred_mat statistics (not future!)
    vals <- pred_mat[, lyr_name]
    
    scaled <- terra::app(
      future_sub[[lyr_name]],
      fun = function(x) {
        (x - mean(vals)) / sdN(vals)
      }
    )
    
    beta <- abs(
      coef_tab$BetaCoefficient[coef_tab$Variable == lyr_name]
    )
    
    scaled * beta
  })
  
  future_rescaled <- terra::rast(out_layers)
  names(future_rescaled) <- layer_names
  future_rescaled
  
}

# ------------------------------------------------------------------------------


#' Cluster novel climate areas independently using NbClust
#'
#' Only called when MESS identifies novel cells.  Samples from novel space,
#' lets NbClust pick the optimal number of clusters, then trains a KNN to
#' paint the novel mask.
#'
#' @param future_rescaled SpatRaster — rescaled + coordinate-weighted predictors.
#' @param novel_mask      Logical SpatRaster identifying novel cells.
#' @param n_novel_pts     Number of points to sample.
#' @param next_cluster_id First integer ID to assign (typically `max(existing) + 1`).
#' @param nbclust_args    User overrides for NbClust.
#'
#' @returns List with single element `clusters_raster` — a SpatRaster covering
#'   only the novel mask extent (NA elsewhere).
#' @keywords internal
#' @noRd
cluster_novel_areas <- function(
  future_rescaled,
  novel_mask,
  n_novel_pts,
  next_cluster_id,
  nbclust_args
) {

  message("Clustering novel climate areas...")

  # Sample from novel areas only
  novel_pts <- terra::spatSample(
    terra::mask(future_rescaled, novel_mask),
    size   = n_novel_pts,
    method = "random",
    na.rm  = TRUE,
    xy     = TRUE,
    as.points = FALSE
  )

  clust_data <- novel_pts[, !names(novel_pts) %in% c("x", "y")]
  clust_data <- clust_data[stats::complete.cases(clust_data), ]

  # Guard: too few points to meaningfully cluster
  if (nrow(clust_data) < 10) {
    warning("Too few novel climate points for clustering. Assigning single novel cluster.")
    novel_clusters <- novel_mask
    novel_clusters <- terra::ifel(novel_mask == 1, next_cluster_id, NA)
    return(
      list(
        clusters_raster = novel_clusters
      )
    )
  }

  # remove highly collinear features for nbclust
  clust_data <- prep_for_nbclust(clust_data)

  # Distance matrix & NbClust
  d <- stats::dist(clust_data, method = "euclidean")

  nbclust_defaults <- list(
    data     = clust_data,
    diss     = d,
    distance = NULL,
    min.nc   = 2,
    max.nc   = 20,
    method   = "ward.D2",
    index    = "all"
  )

  nb_result <- suppressMessages(
    do.call(NbClust::NbClust, modifyList(nbclust_defaults, nbclust_args))
  )

  optimal_k <- max(nb_result$Best.partition)
  # Cut the tree
  hc                <- stats::hclust(d, method = "ward.D2")
  cluster_assignments <- stats::cutree(hc, k = optimal_k) + next_cluster_id - 1

  # Train KNN on novel clusters to paint the raster
  train_data      <- novel_pts[, !names(novel_pts) %in% c("x", "y")]
  train_data$ID   <- factor(cluster_assignments)

  novel_knn <- suppressMessages(
    caret::train(
      ID ~ .,
      data       = train_data,
      method     = "knn",
      trControl  = caret::trainControl(method = "cv", number = 5),
      metric     = "Accuracy"
    )
  )

  # Predict onto masked raster only
  novel_clusters <- terra::predict(
    terra::mask(future_rescaled, novel_mask),
    model = novel_knn,
    na.rm = TRUE
  )

  list(clusters_raster = novel_clusters)
}


# ------------------------------------------------------------------------------


#' Analyze hierarchical relationships between existing and novel clusters
#'
#' Samples an equal number of points from every cluster (existing and novel),
#' builds a single hierarchical tree on that combined sample, then extracts
#' nearest-neighbour relationships via `cluster::silhouette`.
#'
#' @param clusters_raster      Combined SpatRaster (existing + novel IDs).
#' @param future_rescaled      Rescaled predictor SpatRaster (for extraction).
#' @param existing_ids         Integer vector of existing cluster IDs.
#' @param novel_ids            Integer vector of novel cluster IDs.
#' @param n_sample_per_cluster Points to draw from each cluster.
#'
#' @returns Tibble with columns `novel_cluster_id`, `nearest_existing_id`,
#'   `avg_silhouette_width`.
#' @keywords internal
#' @noRd
analyze_cluster_relationships <- function(
  clusters_raster,
  future_rescaled,
  existing_ids,
  novel_ids,
  n_sample_per_cluster
) {

  if (length(novel_ids) == 0) {
    return(tibble::tibble(
      novel_cluster_id    = integer(),
      nearest_existing_id = integer(),
      avg_silhouette_width = numeric()
    ))
  }

  message("Analyzing relationships between existing and novel clusters...")

  # ---- Convert cluster raster to categorical ----
  all_ids <- c(existing_ids, novel_ids)
  
  # Create a categorical version with explicit levels
  clusters_cat <- terra::as.factor(clusters_raster)

  # Stack with predictors
  stack_for_sampling <- c(future_rescaled, clusters_cat)
  names(stack_for_sampling)[terra::nlyr(stack_for_sampling)] <- "cluster_id"

  # Sample (almost) equally from each cluster
  pts <- terra::spatSample(
    clusters_raster,
    size   = n_sample_per_cluster,   
    method = "stratified",
    na.rm  = TRUE,
    xy = TRUE
  ) |>
    as.data.frame() |>
    sf::st_as_sf(coords = c('x', 'y'), crs = terra::crs(clusters_raster))

  clust_data <- terra::extract(future_rescaled, pts, ID = FALSE)
  clust_data <- clust_data[stats::complete.cases(clust_data), ]
  
  # ---- Distance matrix ----
  cluster_ids <- as.integer(pts[['class']])
  
  if (length(cluster_ids) < 2 || length(unique(cluster_ids)) < 2) {
    warning("Insufficient data for silhouette analysis")
    return(tibble::tibble(
      novel_cluster_id     = integer(),
      nearest_existing_id  = integer(),
      avg_silhouette_width = numeric()
    ))
  }
  
  d <- stats::dist(clust_data, method = "euclidean")

  # ---- Silhouette ----
  sil <- cluster::silhouette(cluster_ids, d)

  if (!is.matrix(sil)) {
    warning("Silhouette returned non-matrix result")
    return(tibble::tibble(
    novel_cluster_id     = integer(),
    nearest_existing_id  = integer(),
    avg_silhouette_width = numeric()
    ))
  }

  # Aggregate: for each novel cluster, most common existing neighbor
  data.frame(
    cluster  = sil[, "cluster"],
    neighbor = sil[, "neighbor"],
    sil_width = sil[, "sil_width"]
  ) |>
    dplyr::filter(cluster %in% novel_ids) |>
    dplyr::group_by(novel_cluster_id = cluster) |>
    dplyr::summarise(
      nearest_existing_id = {
        existing_neighbors <- neighbor[neighbor %in% existing_ids]
        if (length(existing_neighbors) > 0) {
          as.integer(names(sort(table(existing_neighbors), decreasing = TRUE)[1]))
        } else {
          NA_integer_
        }
      },
      avg_silhouette_width = mean(sil_width),
      .groups = "drop"
    )
}


# ------------------------------------------------------------------------------

#' Calculate per-cluster area and centroid-shift changes
#'
#' Compares current and future cluster polygons.  Handles clusters that
#' persist, disappear, or are newly gained.
#'
#' @param current_sf  SF object of current clusters (must have `$ID`).
#' @param future_sf   SF object of future clusters (must have `$ID`).
#' @param planar_proj Planar projection for area / distance calculations.
#'
#' @returns Tibble with columns `cluster_id`, `current_area_km2`,
#'   `future_area_km2`, `area_change_pct`, `centroid_shift_km`.
#' @keywords internal
#' @noRd
calculate_changes <- function(current_sf, future_sf, planar_proj) {

  current_p <- sf::st_transform(current_sf, planar_proj)
  future_p  <- sf::st_transform(future_sf,  planar_proj)

  ## calculate total areas of each class
  current_p$area <- as.numeric(sf::st_area(current_p)) / 1e6   # km²
  future_p$area  <- as.numeric(sf::st_area(future_p))  / 1e6

  ## find geometric centroid
  sf::st_agr(current_p) <- 'constant'
  sf::st_agr(future_p) <- 'constant'

  current_cents <- sf::st_point_on_surface(current_p)
  future_cents  <- sf::st_point_on_surface(future_p)

  # ---------- clusters present in both ----------
  common_ids <- intersect(current_sf$ID, future_sf$ID)

  persists <- tibble::tibble(
    cluster_id       = common_ids,
    current_area_km2 = current_p$area[match(common_ids, current_sf$ID)],
    future_area_km2  = future_p$area[match(common_ids, future_sf$ID)],
  ) |>
    dplyr::mutate(
      area_change_pct  = 100 * (future_area_km2 - current_area_km2) / current_area_km2,
      centroid_shift_km = as.numeric(
        sf::st_distance(
          current_cents[match(common_ids, current_sf$ID), ],
          future_cents[match(common_ids, future_sf$ID), ],
          by_element = TRUE
        )
      ) / 1000
    )

  # ---------- lost clusters ----------
  lost_ids <- setdiff(current_sf$ID, future_sf$ID)
  lost <- if (length(lost_ids) > 0) {
    tibble::tibble(
      cluster_id       = lost_ids,
      current_area_km2 = current_p$area[match(lost_ids, current_sf$ID)],
      future_area_km2  = 0,
      area_change_pct  = -100,
      centroid_shift_km = NA_real_
    )
  }

  # ---------- gained clusters ----------
  gained_ids <- setdiff(future_sf$ID, current_sf$ID)
  gained <- if (length(gained_ids) > 0) {
    tibble::tibble(
      cluster_id       = gained_ids,
      current_area_km2 = 0,
      future_area_km2  = future_p$area[match(gained_ids, future_sf$ID)],
      area_change_pct  = NA_real_,
      centroid_shift_km = NA_real_
    )
  }

  # ------- directional movement of point-on-surfaces ------- #

  dplyr::bind_rows(persists, lost, gained) |>
    dplyr::arrange(cluster_id) 
}