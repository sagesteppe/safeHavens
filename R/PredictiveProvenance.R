# ==============================================================================
# projectClusters workflow
# ==============================================================================
# elasticSDM_noPCNM()                  - SDM without spatial autocorrelation
# rescaleFuture()                      - apply current betas to future climate
# projectClusters()                    - main entry point
# cluster_novel_areas()                - NbClust on novel climate only
# analyze_cluster_relationships()      - sample both, force tree, silhouette
# calculate_changes()                  - area / centroid deltas
# ==============================================================================


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
#' @param current_model glmnet model from `elasticSDM_noPCNM()$Model`.
#' @param current_clusters Output list from `EnvironmentalBasedSample()`.
#'   Must contain at minimum `$Geometry` (SF polygons with `$ID`),
#'   `$fit.knn` (trained caret KNN), and access to the training points
#'   used for MESS reference extraction.
#' @param future_predictors SpatRaster of future climate.  Layer names must
#'   match those retained in `current_model`.
#' @param current_predictors SpatRaster of current climate (used for MESS
#'   reference and for `rescaleFuture` standardisation).
#' @param planar_proj EPSG code or proj4 string for a planar projection in
#'   metres (same as used in the current analysis).
#' @param coord_wt Coordinate weighting passed to `add_weighted_coordinates`.
#'   Default `2.5`.
#' @param mess_threshold MESS values below this are treated as novel climate.
#'   Default `0`.
#' @param cluster_novel If `TRUE` and novel cells exist, cluster them
#'   independently with NbClust.  If `FALSE` novel cells are left as `NA`.
#'   Default `TRUE`.
#' @param n_novel_pts Number of points to sample from novel areas for
#'   clustering.  Default `500`.
#' @param n_sample_per_cluster Number of points to sample from each cluster
#'   (existing *and* novel) for the relationship tree.  Default `50`.
#' @param nbclust_args Named list of arguments forwarded to
#'   `NbClust::NbClust`.  Sensible defaults are set internally
#'   (`min.nc = 2`, `max.nc = 10`, `method = "ward.D2"`, `index = "all"`).
#'
#' @returns A named list:
#'   \item{clusters_raster}{`SpatRaster` — cluster IDs across the full extent.}
#'   \item{clusters_sf}{`sf` polygons with column `ID`.}
#'   \item{novel_mask}{`SpatRaster` — logical, `TRUE` where MESS < threshold.}
#'   \item{mess}{`SpatRaster` — raw MESS scores (minimum across all variables).}
#'   \item{changes}{`tibble` — per-cluster area and centroid-shift metrics.}
#'   \item{novel_similarity}{`tibble` — nearest existing cluster and mean
#'     silhouette width for each novel cluster.  Zero-row tibble if none.}
#'
#' @export
projectClusters <- function(
  current_model,
  current_clusters,
  future_predictors,
  current_predictors,
  planar_proj,
  coord_wt            = 2.5,
  mess_threshold      = 0,
  cluster_novel       = TRUE,
  n_novel_pts         = 500,
  n_sample_per_cluster = 50,
  nbclust_args        = list()
) {

  # ---- 1. Rescale future predictors with current betas ----
  future_rescaled <- rescaleFuture(
    model             = current_model,
    future_predictors = future_predictors,
    current_predictors = current_predictors
  )

  # not sure rescaling is working... 


  # ---- 2. MESS — identify novel climate cells ----
  vars <- names(future_rescaled)

  ref_matrix <- current_clusters$TrainData |> 
    dplyr::select(-dplyr::any_of(c('x', 'y', 'ID')))
  ref_matrix <- ref_matrix[stats::complete.cases(ref_matrix), ]

  mess_result  <- terra::rast(  
    dismo::mess(
      x = raster::stack(future_predictors[[vars]]), 
      v = ref_matrix, 
      full = FALSE
    )
  )
  mess_overall <- mess_result[[terra::nlyr(mess_result)]]   # minimum across variables
  novel_mask   <- mess_overall < mess_threshold

  return(mess_result)
  # ---- 3. Add weighted coordinates ----
  future_rescaled <- add_weighted_coordinates(future_rescaled, coord_wt)

  # ---- 4. Predict known-climate cells with existing KNN ----
  known_clusters <- terra::predict(
    future_rescaled,
    model = current_clusters$fit.knn,
    na.rm = TRUE
  )

  # Mask out novel areas — will be filled in step 5 if requested
  known_clusters <- terra::mask(known_clusters, novel_mask, inverse = TRUE)

  # ---- 5. Cluster novel areas independently (two-tree approach) ----
  n_novel_cells <- terra::global(novel_mask, "sum", na.rm = TRUE)[[1]]

  if (cluster_novel && n_novel_cells > 0) {

    novel_result <- cluster_novel_areas(
      future_rescaled  = future_rescaled,
      novel_mask       = novel_mask,
      n_novel_pts      = n_novel_pts,
      next_cluster_id  = max(current_clusters$Geometry$ID) + 1,
      nbclust_args     = nbclust_args
    )

    novel_clusters <- novel_result$clusters_raster

    # Merge: known cells keep their IDs, novel cells get new IDs
    final_clusters <- known_clusters
    final_clusters[!is.na(novel_clusters)] <- novel_clusters[!is.na(novel_clusters)]

    # ---- 6. Relationship analysis (combined tree + silhouette) ----
    novel_similarity <- analyze_cluster_relationships(
      clusters_raster      = final_clusters,
      future_rescaled      = future_rescaled,
      existing_ids         = unique(current_clusters$Geometry$ID),
      novel_ids            = unique(terra::values(novel_clusters, na.rm = TRUE)),
      n_sample_per_cluster = n_sample_per_cluster
    )

  } else {
    final_clusters   <- known_clusters
    novel_similarity <- tibble::tibble(
      novel_cluster_id    = integer(),
      nearest_existing_id = integer(),
      avg_silhouette_width = numeric()
    )
  }

  # ---- 6b. Remove minor cell flecks / noise ----- #
  # Preserve original geometry
  template <- terra::rast(final_clusters)

  agg <- terra::aggregate(
    final_clusters,
    fact = 2,
    fun  = modal,
    na.rm = TRUE
  )

  final_clean <- terra::resample(
    agg,
    template,
    method = "near"
  )
  
  # ---- 7. Polygonise ----
  clusters_sf <- terra::as.polygons(final_clusters) |>
    sf::st_as_sf() |>
    sf::st_make_valid()
  names(clusters_sf)[1] <- "ID"

  # ---- 8. Change metrics ----
  changes <- calculate_changes(
    current_sf = current_clusters$Geometry,
    future_sf  = clusters_sf,
    planar_proj = planar_proj
  )

  # ---- 9. Return ----
  list(
    clusters_raster  = final_clusters,
    clusters_sf      = clusters_sf,
    novel_mask       = novel_mask,
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
#'
#' @returns SpatRaster with rescaled future predictors (one layer per
#'   non-zero coefficient).
#' @export
rescaleFuture <- function(model, future_predictors, current_predictors) {

  # Extract non-zero coefficients (glmnet shrinks some to zero)
  coefs <- as.matrix(stats::coef(model))
  coefs <- coefs[coefs[, 1] != 0, , drop = FALSE]

  # Variable names — drop intercept row
  vars <- rownames(coefs)[-1]

  # Subset both stacks to the same variables
  future_sub  <- future_predictors[[vars]]
  current_sub <- current_predictors[[vars]]

  # Standardisation params from current climate
  current_vals <- terra::as.data.frame(current_sub, na.rm = TRUE)
  means <- colMeans(current_vals)
  sds   <- apply(current_vals, 2, stats::sd)
  sds[sds == 0] <- 1

  # Standardise future using current params, then weight by betas
  # Both operations are vectorised over layers
  for (i in seq_along(vars)) {
    future_sub[[i]] <- (future_sub[[i]] - means[i]) / sds[i] * coefs[vars[i], 1]
  }

  future_sub
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
    novel_clusters[novel_mask == TRUE] <- next_cluster_id
    return(list(clusters_raster = novel_clusters))
  }

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

  all_ids <- c(existing_ids, novel_ids)

  # Sample (almost) equally from each cluster
  sampled <- terra::spatSample(
    c(future_rescaled, clusters_raster),
    size   = n_sample_per_cluster,   # per stratum, each=TRUE is default
    method = "stratified",
    strata = clusters_raster,
    na.rm  = TRUE
  )

  all_samples <- dplyr::bind_rows(sampled)
  all_samples <- all_samples[stats::complete.cases(all_samples), ]

  # Distance matrix on combined sample
  clust_data  <- all_samples[, !names(all_samples) %in% "cluster_id"]
  cluster_ids <- all_samples$cluster_id
  d           <- stats::dist(clust_data, method = "euclidean")

  # Silhouette gives us neighbor cluster for each point
  sil <- cluster::silhouette(cluster_ids, d)

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
 # lwgeom::st_geod_azimuth()st_centroid

  dplyr::bind_rows(persists, lost, gained) |>
    dplyr::arrange(cluster_id)
}