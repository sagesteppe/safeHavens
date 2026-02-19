#' Fit a Bayesian spatial GLMM as a species distribution model
#'
#' @description
#' A Bayesian alternative to [elasticSDM()] using `brms`. Fits a binomial GLMM
#' with a Gaussian process spatial random effect over the full pairwise distance
#' matrix (via a kernel covariance, no mesh approximation). All environmental
#' predictors are retained; shrinkage toward zero is handled by priors rather
#' than variable selection.
#'
#' Users should either (a) remove highly correlated predictors beforehand, or
#' (b) pass PCA axes as predictors (see `pca_predictors` argument), which is
#' the recommended workflow.
#'
#' The function returns an object structurally compatible with the `safeHavens`
#' downstream workflow: the same named list slots are populated where
#' meaningful, with Bayesian-specific additions (posterior draws, LOO-CV).
#'
#' @param x An `sf` object of occurrence points. Must have an `occurrence`
#'   column coded 0/1, OR the function will assume all rows are presences and
#'   generate pseudo-absences (matching [elasticSDM()] behaviour).
#' @param predictors A `terra` SpatRaster stack of environmental predictors.
#' @param planar_projection Numeric EPSG code or proj4 string for a planar
#'   (metre-unit) CRS. Used for spatial thinning and to derive GP coordinates.
#'   See [elasticSDM()] for guidance.
#' @param quantile_v Numeric thinning quantile passed to [spThin::thin()].
#'   Defaults to `0.025`. Set to `0.001` for essentially no thinning.
#' @param prior_type Character. Prior family placed on the fixed-effect
#'   environmental coefficients. One of:
#'   \describe{
#'     \item{`"horseshoe"`}{Regularised horseshoe (Piironen & Vehtari 2017).
#'       Best when many predictors may be irrelevant. Requires the expected
#'       number of non-zero coefficients via `p0` (see `...`).}
#'     \item{`"normal"`}{`Normal(0, prior_scale)`. Simple and fast; reasonable
#'       when predictors are pre-selected or are PCA axes.}
#'     \item{`"student"`}{Student-t(3, 0, prior_scale). Heavier tails than
#'       Normal; a middle ground.}
#'   }
#'   Defaults to `"horseshoe"`.
#' @param prior_scale Numeric. Scale parameter for `"normal"` and `"student"`
#'   priors. Ignored for `"horseshoe"`. Defaults to `1` (weakly informative on
#'   the log-odds scale when predictors are standardised).
#' @param pca_predictors Logical. If `TRUE`, replace the raw predictor stack
#'   with its principal components before fitting. The number of axes retained
#'   is controlled by `pca_axes`. Defaults to `FALSE`.
#' @param pca_axes Integer. Number of PCA axes to retain when
#'   `pca_predictors = TRUE`. Defaults to `5`.
#' @param gp_scale_prior A `brms` prior object for the GP length-scale
#'   parameter. Defaults to `prior(inv_gamma(3.5, 0.5), class = lscale)`, which
#'   places most mass on moderate spatial ranges and is weakly informative for
#'   ecological data. Override if your species has a known range.
#' @param chains Integer. Number of MCMC chains. Defaults to `4`.
#' @param iter Integer. Total iterations per chain (including warmup). Defaults
#'   to `2000`.
#' @param warmup Integer. Warmup iterations per chain. Defaults to `1000`.
#' @param cores Integer. Parallel cores. Defaults to `parallel::detectCores()`.
#' @param k Integer. Number of spatial CV folds (CAST `knndm`). Defaults to
#'   `5`.
#' @param seed Integer. Random seed for reproducibility. Defaults to `42`.
#' @param backend Character. Stan backend passed to [brms::brm()]. One of
#'   `"cmdstanr"` (recommended, faster) or `"rstan"`. Defaults to
#'   `"cmdstanr"`.
#' @param ... Additional arguments forwarded to [brms::brm()]. Also accepts
#'   `p0` (integer, expected non-zero predictors) for the horseshoe prior;
#'   defaults to `floor(ncol(pred_matrix) / 2)`.
#' @param fact Numeric, default 2.0.
#'  Factor to multiple the number of occurrence records by to generate the number of background (absence) points. 

#'
#' @returns A named list with the following elements, mirroring [elasticSDM()]:
#' \describe{
#'   \item{`RasterPredictions`}{`SpatRaster` of posterior mean predicted
#'     occurrence probability.}
#'   \item{`RasterPredictions_sd`}{`SpatRaster` of posterior SD (uncertainty).}
#'   \item{`Predictors`}{The predictor raster stack actually used (may be PCA
#'     axes if `pca_predictors = TRUE`).}
#'   \item{`PCNM`}{Always `NULL`; slot retained for workflow compatibility.}
#'   \item{`Model`}{The fitted `brmsfit` object.}
#'   \item{`CVStructure`}{CAST `knndm` object with fold indices.}
#'   \item{`LOO`}{`loo` object from approximate leave-one-out cross-validation
#'     via PSIS-LOO. Replaces `ConfusionMatrix`.}
#'   \item{`ConfusionMatrix`}{A simple threshold-based confusion matrix
#'     (threshold = 0.5 on posterior mean) on the held-out test fold, for
#'     comparability with [elasticSDM()] output.}
#'   \item{`TrainData`}{`sf` object used for model fitting.}
#'   \item{`TestData`}{`sf` object held out for evaluation.}
#'   \item{`PredictMatrix`}{The numeric matrix of predictor values used for
#'     the spatial raster prediction.}
#'   \item{`PCAModel`}{The `prcomp` object if `pca_predictors = TRUE`,
#'     otherwise `NULL`. Retained so users can project new data onto the same
#'     axes.}
#'   \item{`Diagnostics`}{Named list: `Rhat` (max R-hat across all parameters),
#'     `BulkESS` (min bulk ESS), `TailESS` (min tail ESS). Flags convergence
#'     problems automatically.}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' x <- read.csv(file.path(system.file(package = "dismo"), "ex", "bradypus.csv"))
#' x <- sf::st_as_sf(x[, c("lon", "lat")], coords = c("lon", "lat"), crs = 4326)
#'
#' files <- list.files(
#'   path    = file.path(system.file(package = "dismo"), "ex"),
#'   pattern = "grd",
#'   full.names = TRUE
#' )
#' predictors <- terra::rast(files)
#'
#' # Recommended workflow: PCA axes as predictors, horseshoe prior
#' result <- bayesianSDM(
#'   x                 = x,
#'   predictors        = predictors,
#'   planar_projection = 5070,
#'   pca_predictors    = TRUE,
#'   pca_axes          = 5,
#'   prior_type        = "horseshoe",
#'   chains            = 4,
#'   iter              = 2000
#' )
#'
#' terra::plot(result$RasterPredictions)
#' terra::plot(result$RasterPredictions_sd)  # uncertainty surface
#' print(result$Diagnostics)
#' }
#'
#' @seealso [elasticSDM()] for the elastic-net alternative,
#'   [brms::brm()], [CAST::knndm()]
#' @export
bayesianSDM <- function(
    x,
    predictors,
    planar_projection,
    quantile_v        = 0.025,
    prior_type        = c("horseshoe", "normal", "student"),
    prior_scale       = 1,
    pca_predictors    = TRUE,
    pca_axes          = 5,
    gp_scale_prior    = brms::prior(inv_gamma(3.5, 0.5), class = lscale),
    chains            = 4,
    iter              = 2000,
    warmup            = 1000,
    cores             = parallel::detectCores(),
    k                 = 5,
    seed              = 42,
    backend           = "cmdstanr",
    fact              = 2,
    ...
) {
  prior_type <- match.arg(prior_type)
  dots       <- list(...)

  # --- 0. PCA on raster surface -----------------
  ## note doing this first with random points allows for clean eDist sampling in the presence 
  # of autocorrelated values. 
  pca_model <- NULL
  if (pca_predictors) {
    message(sprintf("Running PCA — retaining %d axes.", pca_axes))
    pca_result <- run_pca_on_predictors(x, predictors, pca_axes)
    x          <- pca_result$data          # sf with PC columns
    predictors <- pca_result$raster        # SpatRaster of PC axes
    pred_names <- pca_result$pc_names
    pca_model  <- pca_result$pca_model
  }

  # ── 1. Background points & occurrence column ────────────────────────────────
  if (!"occurrence" %in% names(x)) {
    message("No `occurrence` column found — generating pseudo-absences.")
    x$occurrence <- 1L
    pa <- generate_background_points(predictors, x, fact)
    x  <- dplyr::bind_rows(x, pa)
  }
  x <- dplyr::mutate(x, occurrence = as.integer(as.character(occurrence)))

  # ── 2. Spatial thinning ─────────────────────────────────────────────────────
  x <- thin_occurrence_points(x, quantile_v)

  # ── 3. Extract predictor values ─────────────────────────────────────────────
  x <- extract_predictors_to_points(x, predictors)

  # Drop any rows with NA in extracted predictors
  pred_names <- names(predictors)
  complete   <- stats::complete.cases(sf::st_drop_geometry(x)[, pred_names])
  x          <- x[complete, ]

  # ── 5. Train / test split ────────────────────────────────────────────────────
  index <- unlist(caret::createDataPartition(x$occurrence, p = 0.8))
  train <- x[ index, ]
  test  <- x[-index, ]

  # ── 6. Spatial CV folds (CAST knndm) ────────────────────────────────────────
  cv_folds <- create_spatial_cv_folds(train, predictors, k = k)

  # ── 7. Build GP coordinates (planar, scaled) ─────────────────────────────────
  # brms gp() uses raw coordinate columns; we work in km for numerical stability
  train <- add_planar_coords(train, planar_projection, scale_km = TRUE)
  test  <- add_planar_coords(test,  planar_projection, scale_km = TRUE)

  # ── 8. Assemble model formula ────────────────────────────────────────────────
  env_terms <- paste(pred_names, collapse = " + ")
  bf_formula <- stats::as.formula(
    sprintf(
      "occurrence ~ %s + gp(gp_x, gp_y, scale = FALSE, cov = 'exp_quad')",
      env_terms
    )
  )

  # ── 9. Priors ────────────────────────────────────────────────────────────────
  model_prior <- build_priors(
    prior_type  = prior_type,
    prior_scale = prior_scale,
    pred_names  = pred_names,
    n_obs       = nrow(train),
    p0          = dots[["p0"]],
    gp_prior    = gp_scale_prior
  )

  # ── 10. Fit model ─────────────────────────────────────────────────────────────
  message("Fitting Bayesian spatial GLMM via brms ...")
  fit <- brms::brm(
    formula  = bf_formula,
    data     = sf::st_drop_geometry(train),
    family   = brms::bernoulli(link = "logit"),
    prior    = model_prior,
    chains   = chains,
    iter     = iter,
    warmup   = warmup,
    cores    = cores,
    seed     = seed,
    backend  = backend,
    ...
  )

  # ── 11. Convergence diagnostics ───────────────────────────────────────────────
  diagnostics <- check_convergence(fit)
  if (diagnostics$max_Rhat > 1.05) {
    warning(
      sprintf(
        "Max R-hat = %.3f > 1.05. Chains may not have converged. ",
        diagnostics$max_Rhat
      ),
      "Consider increasing `iter` or inspecting pairs/trace plots."
    )
  }

  # ── 12. LOO cross-validation ──────────────────────────────────────────────────
  message("Computing PSIS-LOO ...")
  loo_result <- brms::loo(fit, moment_match = TRUE)

  # ── 13. Test-set confusion matrix ─────────────────────────────────────────────
  cm <- evaluate_bayes_model(fit, test, pred_names)

  # -- 14. Area of Applicability surface
  message("Computing Area of Applicability (AOA) ...")
  aoa_result <- CAST::aoa(
    train      = sf::st_drop_geometry(train[, pred_names]),
    predictors = predictors[[pred_names]],
    model      = fit
  )

  # ── 14. Spatial raster predictions ───────────────────────────────────────────
  message("Generating raster predictions (posterior mean & SD) ...")
  rast_list <- create_bayes_spatial_predictions(fit, predictors, pred_names,
                                                planar_projection)

  list(
    RasterPredictions    = rast_list$mean,
    RasterPredictions_sd = rast_list$sd,
    Predictors           = predictors,
    PCNM                 = NULL,
    Model                = fit,
    CVStructure          = cv_folds,
    LOO                  = loo_result,
    ConfusionMatrix      = cm,
    TrainData            = train,
    TestData             = test,
    PredictMatrix        = rast_list$pred_matrix,
    AOA        = aoa_result$AOA,        # Binary raster: 1 = inside AOA
    DI         = aoa_result$DI,         # Dissimilarity index raster
    AOA_Diagnostics = list(
      max_Rhat    = ...,
      AOA_coverage = sum(terra::values(aoa_result$AOA), na.rm = TRUE) / 
                     sum(!is.na(terra::values(aoa_result$AOA)))
    ), 
    PCAModel             = pca_model,
    Diagnostics          = diagnostics
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

#' Add scaled planar coordinates to an sf object
#'
#' Projects to `planar_projection`, extracts XY coordinates in kilometres,
#' and appends as columns `gp_x` and `gp_y`. Scaling to km avoids the
#' numerical issues that arise when GP length-scales must span millions of
#' metres.
#'
#' @param x sf object
#' @param planar_projection EPSG or proj4 string
#' @param scale_km Logical; divide metres by 1000? Defaults to TRUE.
#' @return sf object with `gp_x` and `gp_y` columns
#' @keywords internal
#' @noRd
add_planar_coords <- function(x, planar_projection, scale_km = TRUE) {
  coords <- sf::st_coordinates(sf::st_transform(x, planar_projection))
  divisor <- if (scale_km) 1000 else 1
  x$gp_x <- coords[, "X"] / divisor
  x$gp_y <- coords[, "Y"] / divisor
  x
}


#' Run PCA on raster predictors and project both points and raster
#'
#' Fits PCA on the raster cell values (a random sample of up to 2500 cells
#' to keep memory tractable), then projects both the point data and the full
#' raster onto the retained axes.
#'
#' @param x sf object with predictor columns already extracted
#' @param predictors SpatRaster
#' @param pca_axes Integer, number of axes to retain
#' @return List: `data` (sf with PC columns), `raster` (SpatRaster of PCs),
#'   `pc_names` (character), `pca_model` (prcomp)
#' @keywords internal
#' @noRd
run_pca_on_predictors <- function(x, predictors, pca_axes) {
  # Fit PCA on raster values (sample for speed)
  pred_names <- names(predictors)
  rast_vals <- terra::spatSample(
    predictors,
    size    = min(2500, terra::ncell(predictors)),
    method  = "random",
    na.rm   = TRUE
  )
  pca_model <- stats::prcomp(rast_vals, center = TRUE, scale. = TRUE)

  # Number of axes cannot exceed number of variables
  pca_axes <- min(pca_axes, length(pred_names))
  pc_names <- paste0("PC", seq_len(pca_axes))

  # Project point data
  point_mat  <- as.matrix(sf::st_drop_geometry(x)[, pred_names])
  point_pcs  <- stats::predict(pca_model, newdata = point_mat)[, seq_len(pca_axes), drop = FALSE]
  x_pc       <- x[, setdiff(names(x), pred_names)]  # drop raw predictor cols
  for (nm in pc_names) x_pc[[nm]] <- point_pcs[, nm]

  # Project raster — process layer by layer to avoid memory blow-up
  # We use terra::predict which calls stats::predict on a per-cell basis
  pc_raster <- terra::predict(
    predictors[[pred_names]],
    model = pca_model,
    index = seq_len(pca_axes)
  )
  names(pc_raster) <- pc_names

  list(
    data      = x_pc,
    raster    = pc_raster,
    pc_names  = pc_names,
    pca_model = pca_model
  )
}


#' Build brms prior set
#'
#' @param prior_type One of "horseshoe", "normal", "student"
#' @param prior_scale Numeric scale for normal/student
#' @param pred_names Character vector of predictor names
#' @param n_obs Number of training observations
#' @param p0 Expected non-zero coefficients (horseshoe only)
#' @param gp_prior brms prior for GP length-scale
#' @return brms `brmsprior` object
#' @keywords internal
#' @noRd
build_priors <- function(prior_type, prior_scale, pred_names, n_obs, p0, gp_prior) {
  p <- length(pred_names)

  env_prior <- switch(prior_type,
    horseshoe = {
      # Regularised horseshoe: Piironen & Vehtari (2017) default parameterisation
      # p0: expected number of non-zero slopes; if unset, use half of predictors
      p0_val <- if (is.null(p0)) floor(p / 2) else as.integer(p0)
      brms::prior(horseshoe(df = 1, par_ratio = NULL), class = b)
      # Note: brms handles the scale via global_scale = p0 / (p - p0) / sqrt(n)
      # exposed via set_prior() string form for full control:
      brms::set_prior(
        sprintf(
          "horseshoe(df = 1, scale_global = %f, df_global = 1, scale_slab = 2, df_slab = 4)",
          (p0_val / (p - p0_val + 0.5)) / sqrt(n_obs)
        ),
        class = "b"
      )
    },
    normal  = brms::set_prior(sprintf("normal(0, %f)", prior_scale),  class = "b"),
    student = brms::set_prior(sprintf("student_t(3, 0, %f)", prior_scale), class = "b")
  )

  intercept_prior <- brms::set_prior("normal(0, 2.5)", class = "Intercept")

  # GP marginal SD: half-normal weakly informative
  sdgp_prior <- brms::set_prior("normal(0, 1)", class = "sdgp")

  c(env_prior, intercept_prior, sdgp_prior, gp_prior)
}


#' Compute MCMC convergence diagnostics
#'
#' @param fit brmsfit object
#' @return Named list: max_Rhat, min_BulkESS, min_TailESS
#' @keywords internal
#' @noRd
check_convergence <- function(fit) {
  summ <- brms::posterior_summary(fit)
  rhats <- brms::rhat(fit)
  ess_b <- brms::neff_ratio(fit)  # relative ESS

  # Raw ESS from posterior package
  draws   <- brms::as_draws_df(fit)
  ess_raw <- posterior::summarise_draws(
    draws,
    bulk_ess = posterior::ess_bulk,
    tail_ess = posterior::ess_tail
  )

  list(
    max_Rhat    = max(rhats, na.rm = TRUE),
    min_BulkESS = min(ess_raw$bulk_ess, na.rm = TRUE),
    min_TailESS = min(ess_raw$tail_ess, na.rm = TRUE),
    n_divergent = sum(brms::nuts_params(fit)$Value[
      brms::nuts_params(fit)$Parameter == "divergent__"
    ])
  )
}


#' Evaluate model on held-out test data and return a simple confusion matrix
#'
#' Posterior mean predicted probability is thresholded at 0.5.
#'
#' @param fit brmsfit
#' @param test_data sf object with `occurrence` and `gp_x`/`gp_y` columns
#' @param pred_names Character vector of predictor column names
#' @return caret confusionMatrix object
#' @keywords internal
#' @noRd
evaluate_bayes_model <- function(fit, test_data, pred_names) {
  pred_df <- sf::st_drop_geometry(test_data)

  # epred_rvars returns posterior distribution; rowMeans of draws = posterior mean
  epred <- brms::posterior_epred(
    fit,
    newdata   = pred_df,
    allow_new_levels = TRUE  # needed because GP is a new prediction location
  )
  # epred is draws x observations matrix
  prob_mean <- colMeans(epred)
  predicted <- factor(as.integer(prob_mean >= 0.5), levels = c(0, 1))
  observed  <- factor(as.integer(pred_df$occurrence),  levels = c(0, 1))

  caret::confusionMatrix(
    data      = predicted,
    reference = observed,
    positive  = "1"
  )
}


#' Generate posterior mean and SD raster predictions
#'
#' For each raster cell, obtains the full posterior predictive distribution
#' of the occurrence probability, then summarises to mean and SD rasters.
#' Processes the raster in row chunks to avoid exhausting memory.
#'
#' @param fit brmsfit
#' @param predictors SpatRaster (possibly PCA)
#' @param pred_names Character vector of predictor column names
#' @param planar_projection EPSG or proj4 string
#' @param chunk_size Integer, cells per chunk. Defaults to 10 000.
#' @return List: `mean` SpatRaster, `sd` SpatRaster, `pred_matrix` data.frame
#' @keywords internal
#' @noRd
create_bayes_spatial_predictions <- function(
    fit,
    predictors,
    pred_names,
    planar_projection,
    chunk_size = 10000
) {
  # Pull predictor raster values + xy coordinates
  pred_vals <- as.data.frame(predictors[[pred_names]], xy = TRUE, na.rm = FALSE)
  has_data  <- stats::complete.cases(pred_vals[, pred_names])

  # Planar coordinates in km for the GP term
  coords_geo <- terra::project(
    terra::vect(pred_vals[, c("x", "y")], crs = terra::crs(predictors),
                geom = c("x", "y")),
    planar_projection
  )
  coords_km <- terra::crds(coords_geo) / 1000
  pred_vals$gp_x <- coords_km[, 1]
  pred_vals$gp_y <- coords_km[, 2]

  # Initialise output vectors
  mean_vec <- rep(NA_real_, nrow(pred_vals))
  sd_vec   <- rep(NA_real_, nrow(pred_vals))

  # Chunk prediction to avoid memory overload
  valid_idx <- which(has_data)
  chunks    <- split(valid_idx, ceiling(seq_along(valid_idx) / chunk_size))

  message(sprintf(
    "Predicting over %d valid cells in %d chunks ...",
    length(valid_idx), length(chunks)
  ))

  for (i in seq_along(chunks)) {
    if (i %% 10 == 0 || i == length(chunks)) {
      message(sprintf("  Chunk %d / %d", i, length(chunks)))
    }
    idx   <- chunks[[i]]
    chunk <- pred_vals[idx, c(pred_names, "gp_x", "gp_y")]

    epred <- brms::posterior_epred(
      fit,
      newdata          = chunk,
      allow_new_levels = TRUE
    )
    mean_vec[idx] <- colMeans(epred)
    sd_vec[idx]   <- apply(epred, 2, stats::sd)
  }

  # Rebuild rasters
  template           <- predictors[[1]]
  rast_mean          <- terra::rast(template)
  rast_sd            <- terra::rast(template)
  terra::values(rast_mean) <- mean_vec
  terra::values(rast_sd)   <- sd_vec
  names(rast_mean)   <- "occurrence_prob_mean"
  names(rast_sd)     <- "occurrence_prob_sd"

  list(
    mean         = rast_mean,
    sd           = rast_sd,
    pred_matrix  = pred_vals[, c(pred_names, "gp_x", "gp_y")]
  )
}
