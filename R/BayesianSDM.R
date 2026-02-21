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
#' @param resample Boolean, Defaults to FALSE. Used to place 15% of the requested points in areas undersampled by sdm::background functions. 
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
#'   parameter, or `NULL` to use brms default. The default `NULL` uses
#'   `inv_gamma(3, 1)`, which is weakly informative and appropriate for most
#'   ecological datasets. Override only if you have strong prior knowledge about
#'   the spatial range of your species (e.g., `prior(inv_gamma(5, 2), class = lscale, coef = "")`
#'   for tighter spatial autocorrelation).
#' #' @param feature_selection Character. Variable selection method to apply before
#'   Bayesian model fitting. One of:
#'   \describe{
#'     \item{`"ffs"`}{Forward feature selection via CAST::ffs(). Uses spatial CV
#'       folds to select variables. Fast and spatially-aware. Default.}
#'     \item{`"none"`}{No feature selection; use all predictors (or all PCA axes).
#'       Relies on horseshoe prior for shrinkage.}
#'   }
#' @param min_ffs_var Integer. Minium number of ffs vars to start with. 
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
#'  Factor to multiple the number of occurrence records by to generate the number of background (absence) points. #'
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
    gp_scale_prior    = NULL,
    resample          = FALSE,
    feature_selection = c("ffs", "none"),
    min_ffs_var       = 5, 
    chains            = 4,
    iter              = 2000,
    warmup            = 1000,
    cores             = 4,
    k                 = 5,
    seed              = 42,
    backend           = "cmdstanr",
    fact              = 3, 
    ...
) {
  prior_type <- match.arg(prior_type)
  dots       <- list(...)

  # ── 0. Optional PCA on raster surface ──────────────────────────────────────
  # Run PCA BEFORE background generation so eDist method works on orthogonal axes
  # rather than potentially collinear raw predictors (avoids singular covariance)
  pca_model <- NULL
  if (pca_predictors) {
    pca_result  <- run_pca_on_raster(predictors, pca_axes)
    predictors  <- pca_result$raster        # SpatRaster of PC axes
    pred_names  <- pca_result$pc_names
    pca_model   <- pca_result$pca_model
  } else {
    pred_names <- names(predictors)
  }

  if(min_ffs_var < terra::nlyr(predictors)){min_ffs_var = 2}

  # ── 1. Background points & occurrence column ────────────────────────────────
  if (!"occurrence" %in% names(x)) {
    x$occurrence <- 1L
    pa <- generate_background_points(predictors, x, fact, resample)
    x  <- dplyr::bind_rows(x, pa)
  }
  x <- dplyr::mutate(x, occurrence = as.integer(as.character(occurrence)))

  # ── 2. Spatial thinning ─────────────────────────────────────────────────────
  x <- thin_occurrence_points(x, quantile_v)

  # ── 3. Extract predictor values ─────────────────────────────────────────────
  x <- extract_predictors_to_points(x, predictors)

  # Drop any rows with NA in extracted predictors
  complete   <- stats::complete.cases(sf::st_drop_geometry(x)[, pred_names])
  x          <- x[complete, ]

  # ── 4. Train / test split ────────────────────────────────────────────────────
  index <- unlist(caret::createDataPartition(factor(x$occurrence), p = 0.8))
  train <- x[ index, ]
  test  <- x[-index, ]

  # ── 5. Spatial CV folds (CAST knndm) ────────────────────────────────────────
  cv_folds <- create_spatial_cv_folds(train, predictors, k = k)

    # ── 6. Optional feature selection ───────────────────────────────────────────
  feature_selection <- match.arg(feature_selection)
  if (feature_selection != "none") {
    message(sprintf("Running %s feature selection ...", 
                    toupper(feature_selection)))
    selected_vars <- perform_feature_selection_bayes(
      train          = train,
      pred_names     = pred_names,
      cv_folds       = cv_folds,
      min_ffs_var    = min_ffs_var,
      method         = feature_selection
    )
    
    if (length(selected_vars) == 0) {
      warning("Feature selection returned 0 variables. Using all predictors.")
    } else {
      message(sprintf("  Selected %d of %d variables: %s",
                      length(selected_vars), length(pred_names),
                      paste(selected_vars, collapse = ", ")))
      pred_names <- selected_vars
      predictors <- predictors[[pred_names]]  # Subset raster
    }
  }

  # ── 6. Build GP coordinates (planar, scaled) ─────────────────────────────────
  # brms gp() uses raw coordinate columns; we work in km for numerical stability
  train <- add_planar_coords(train, planar_projection, scale_km = TRUE)
  test  <- add_planar_coords(test,  planar_projection, scale_km = TRUE)

  # ── 7. Assemble model formula ────────────────────────────────────────────────
  env_terms <- paste(pred_names, collapse = " + ")
  bf_formula <- stats::as.formula(
    sprintf(
      "occurrence ~ %s + s(gp_x, gp_y, bs = 'gp', k = 50)",
      env_terms
    )
  )

  # ── 8. Priors ────────────────────────────────────────────────────────────────
  model_prior <- build_priors(
    prior_type  = prior_type,
    prior_scale = prior_scale,
    pred_names  = pred_names,
    n_obs       = nrow(train),
    p0          = dots[["p0"]],
    gp_prior    = gp_scale_prior
  )

  # ── 9. Fit model ─────────────────────────────────────────────────────────────  
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
    save_pars = brms::save_pars(all = TRUE),
    control = list(adapt_delta = 0.99),
    ...
  )

  # ── 10. Convergence diagnostics ───────────────────────────────────────────────
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

  # ── 11. LOO cross-validation ──────────────────────────────────────────────────
  message("Computing PSIS-LOO ...")
  loo_result <- brms::loo(fit, moment_match = TRUE, reloo = TRUE)

  # ── 12. Test-set confusion matrix ─────────────────────────────────────────────
  cm <- evaluate_bayes_model(fit, test, pred_names)

  # ── 13. Spatial raster predictions ───────────────────────────────────────────
  message("Generating raster predictions (posterior mean & SD) ...")
  rast_list <- create_bayes_spatial_predictions(fit, predictors, 
    planar_projection = paste0('EPSG:', planar_projection), pred_names)

  # -- 14. Area of Applicability surface
 # message("Computing Area of Applicability (AOA) ...")
 # aoa_result <- CAST::aoa(
 #   newdata = predictors,
 #   model      = fit,
 #   train      = sf::st_drop_geometry(train[, pred_names])
 # )

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
  #  AOA                  = aoa_result$AOA,       
  #  DI                   = aoa_result$DI,    
 #   AOA_Diagnostics      = list(
 #     max_Rhat           = ...,
 #     AOA_coverage       = sum(terra::values(aoa_result$AOA), na.rm = TRUE) / 
 #                          sum(!is.na(terra::values(aoa_result$AOA)))
 #   ), 
    PCAModel             = pca_model,
    Diagnostics          = diagnostics
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

#' Perform feature selection before Bayesian model fitting
#'
#' Uses either CAST forward feature selection (spatially-aware)
#' to reduce predictor set before expensive Bayesian fitting.
#'
#' @param train sf object with training data
#' @param pred_names Character vector of predictor names
#' @param cv_folds CAST knndm object with CV fold indices
#' @param method Currently only "ffs"
#' @param min_ffs_var Minimum number of FFS variables to start with. 
#' @return Character vector of selected predictor names
#' @keywords internal
#' @noRd
perform_feature_selection_bayes <- function(train, pred_names, cv_folds, method, min_ffs_var) {
  train_df <- sf::st_drop_geometry(train)
  
  if (method == "ffs") {
    # Forward feature selection using spatial CV folds
    ffs_result <- suppressMessages(
      CAST::ffs(
        predictors = train_df[, pred_names, drop = FALSE],
        response   = factor(train_df$occurrence),
        method     = "glm",
        minVar     =  min_ffs_var, 
        family     = stats::binomial(),
        metric     = "Accuracy",
        trControl  = caret::trainControl(
          method = "cv",
          index  = cv_folds$indx_train,
          indexOut = cv_folds$indx_test,
          savePredictions = FALSE,
          verboseIter = FALSE
        ),
        verbose    = FALSE
      )
    )
    return(ffs_result$selectedvars)
  }
  
  # Fallback
  pred_names
}

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


#' Run PCA on raster predictors only
#'
#' Fits PCA on a random sample of raster cells (up to 50,000 for memory),
#' then projects the full raster onto the retained axes. This happens BEFORE
#' any point-based operations so that background sampling (eDist method) works
#' on orthogonal PC axes rather than potentially collinear raw predictors.
#'
#' @param predictors SpatRaster of raw environmental variables
#' @param pca_axes Integer, number of axes to retain
#' @return List: `raster` (SpatRaster of PCs), `pc_names` (character),
#'   `pca_model` (prcomp object)
#' @keywords internal
#' @noRd
run_pca_on_raster <- function(predictors, pca_axes) {
  pred_names <- names(predictors)
  
  # Fit PCA on raster values (sample for speed)
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

  # Project raster — terra::predict calls stats::predict on a per-cell basis
  pc_raster <- terra::predict(
    predictors,
    model = pca_model,
    index = seq_len(pca_axes)
  )
  names(pc_raster) <- pc_names

  list(
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

  c(env_prior, intercept_prior)
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
#' Uses terra::predict() with a custom wrapper function to efficiently
#' predict across the full raster. Much faster than manual chunking because
#' terra handles memory management and parallelization internally.
#'
#' @param fit brmsfit with s(gp_x, gp_y) spatial smooth (not gp() term)
#' @param predictors SpatRaster (possibly PCA)
#' @param pred_names Character vector of predictor column names
#' @param planar_projection EPSG or proj4 string
#' @return List: `mean` SpatRaster, `sd` SpatRaster, `pred_matrix` data.frame
#' @keywords internal
#' @noRd
create_bayes_spatial_predictions <- function(
    fit,
    predictors,
    pred_names,
    planar_projection
) {
  
  # Build coordinate rasters in planar projection (km scale)
  template <- predictors[[1]]
  coords_lonlat <- terra::as.data.frame(template, xy = TRUE, cells = FALSE)[, c("x", "y")]
  
  # Project to planar km
  coords_vect <- terra::vect(
    coords_lonlat,
    geom = c("x", "y"),
    crs = terra::crs(predictors)
  )
  coords_proj <- terra::project(coords_vect, planar_projection)
  coords_km <- terra::crds(coords_proj) / 1000
  
  # Create coordinate rasters
  gp_x_rast <- terra::rast(template)
  gp_y_rast <- terra::rast(template)
  terra::values(gp_x_rast) <- coords_km[, 1]
  terra::values(gp_y_rast) <- coords_km[, 2]
  names(gp_x_rast) <- "gp_x"
  names(gp_y_rast) <- "gp_y"
  
  # Combine predictors with coordinate rasters
  pred_stack <- c(predictors[[pred_names]], gp_x_rast, gp_y_rast)
  
  # Custom predict function for terra
  # terra passes data as a data.frame with column names matching raster layer names
  predict_mean_sd <- function(model, data, ...) {
    # Subsample posterior for speed 
    epred <- brms::posterior_epred(
      model,
      newdata = data,
      allow_new_levels = TRUE,
      ndraws = 2000
    )
    # Return both mean and sd as a 2-column matrix
    cbind(
      mean = colMeans(epred),
      sd   = apply(epred, 2, stats::sd)
    )
  }
  
  # Let terra handle all the chunking/memory management
  result <- terra::predict(
    pred_stack,
    model = fit,
    fun = predict_mean_sd,
    na.rm = TRUE,
    cores = 1  # brms isn't thread-safe for this, keep sequential
  )
  
  # Split into separate rasters
  rast_mean <- result[[1]]
  rast_sd   <- result[[2]]
  names(rast_mean) <- "occurrence_prob_mean"
  names(rast_sd)   <- "occurrence_prob_sd"
  
  # Build pred_matrix for downstream compatibility
  pred_vals <- as.data.frame(pred_stack, xy = FALSE, na.rm = FALSE)
  
  list(
    mean         = rast_mean,
    sd           = rast_sd,
    pred_matrix  = pred_vals
  )
}
