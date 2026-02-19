#' Rescale a raster stack to reflect posterior beta coefficients from a brms SDM
#'
#' @description
#' Bayesian analogue of [RescaleRasters()] for use with [bayesianSDM()] output.
#' Each predictor raster is centred and scaled using training-data moments, then
#' multiplied by the absolute value of its posterior beta coefficient, so that
#' the resulting layers are on a common scale weighted by ecological importance.
#' The rescaled stack can be passed directly into [EnvironmentalBasedSample()].
#'
#' Because `brms` uses `autoscale = TRUE` by default for regularising priors
#' (horseshoe, normal, student), the posterior betas are already on a
#' standardised scale and are directly comparable across predictors — no
#' additional correction is needed beyond what [RescaleRasters()] applies to
#' glmnet output.
#'
#' @section Posterior summary choices:
#' The `beta_summary` argument controls which posterior summary statistic is
#' used as the point-estimate beta weight:
#' \describe{
#'   \item{`"mean"`}{Posterior mean. Minimum MSE estimator; default and
#'     recommended for most use cases.}
#'   \item{`"median"`}{Posterior median. More robust when horseshoe priors
#'     produce heavy-tailed marginals for near-zero coefficients.}
#'   \item{`"Q2.5"` / `"Q97.5"`}{Credible interval bounds. Useful for
#'     sensitivity analysis (e.g. conservative lower-bound weighting).}
#' }
#'
#' @section Uncertainty layer:
#' When `include_uncertainty = TRUE` an additional raster layer is appended
#' whose values are the *posterior SD of the linear predictor* at each cell,
#' propagated through all environmental betas. This encodes where the model is
#' most uncertain about environmental conditions, and can be used as an
#' optional extra input to [EnvironmentalBasedSample()] to bias cluster
#' sampling toward uncertain regions. Set `uncertainty_wt` to control its
#' influence relative to environmental layers (analogous to `coord_wt`).
#'
#' @param model A `brmsfit` object from [bayesianSDM()].
#' @param predictors A `SpatRaster` stack. Should be the `$Predictors` element
#'   from [bayesianSDM()] output (i.e. already PCA-rotated if `pca_predictors`
#'   was `TRUE`).
#' @param training_data An `sf` object. The `$TrainData` element from
#'   [bayesianSDM()] output. Used to compute per-variable mean and SD for
#'   centring/scaling the raster.
#' @param pred_mat A data frame or matrix. The `$PredictMatrix` element from
#'   [bayesianSDM()] output. Column names must match environmental predictor
#'   names (GP coordinate columns `gp_x`/`gp_y` are ignored automatically).
#' @param beta_summary Character. Which posterior summary to use as the beta
#'   weight. One of `"mean"`, `"median"`, `"Q2.5"`, `"Q97.5"`.
#'   Defaults to `"mean"`.
#' @param include_uncertainty Logical. Whether to append a layer encoding
#'   posterior uncertainty of the linear predictor. Defaults to `FALSE`.
#' @param uncertainty_wt Numeric. Weight applied to the uncertainty layer
#'   relative to the maximum environmental layer range, analogous to `coord_wt`
#'   in [EnvironmentalBasedSample()]. Only used when
#'   `include_uncertainty = TRUE`. Defaults to `1`.
#'
#' @returns A list with three elements:
#' \describe{
#'   \item{`RescaledPredictors`}{`SpatRaster` of rescaled, beta-weighted
#'     predictor layers. Pass this as `pred_rescale` to
#'     [EnvironmentalBasedSample()].}
#'   \item{`BetaCoefficients`}{Data frame of posterior summaries for all fixed
#'     effects (columns: `Variable`, `Estimate` (mean), `Est.Error` (SD),
#'     `Q2.5`, `Q97.5`, `BetaWeight` (the value actually used for scaling)).
#'     Intercept and GP parameters are excluded.}
#'   \item{`UncertaintyLayer`}{`SpatRaster` of propagated posterior SD, or
#'     `NULL` if `include_uncertainty = FALSE`.}
#' }
#'
#' @seealso [RescaleRasters()] for the glmnet equivalent,
#'   [bayesianSDM()], [EnvironmentalBasedSample()]
#' @export
RescaleRasters_bayes <- function(
    model,
    predictors,
    training_data,
    pred_mat,
    beta_summary      = c("mean", "median", "Q2.5", "Q97.5"),
    include_uncertainty = FALSE,
    uncertainty_wt    = 1
) {
  beta_summary <- match.arg(beta_summary)

  # ── 1. Extract fixed-effect posterior summaries ──────────────────────────────
  coef_tab <- extract_posterior_betas(model, beta_summary)

  # ── 2. Identify environmental predictors (drop GP columns) ──────────────────
  env_vars <- coef_tab$Variable
  # pred_mat may contain gp_x, gp_y — ignore those
  env_vars <- env_vars[env_vars %in% colnames(pred_mat)]
  env_vars <- env_vars[!env_vars %in% c("gp_x", "gp_y")]
  coef_tab <- coef_tab[coef_tab$Variable %in% env_vars, , drop = FALSE]

  if (nrow(coef_tab) == 0) {
    stop(
      "No environmental predictor betas found. ",
      "Check that `pred_mat` columns match the model fixed effects."
    )
  }

  # ── 3. Validate raster and pred_mat alignment ────────────────────────────────
  missing_rast <- setdiff(env_vars, names(predictors))
  if (length(missing_rast) > 0) {
    stop(
      "These predictors are in the model but not in `predictors` raster: ",
      paste(missing_rast, collapse = ", ")
    )
  }

  # ── 4. Compute training-data moments for centring/scaling ───────────────────
  # Use pred_mat (same data that entered the model) for consistency
  pm <- as.data.frame(pred_mat)[, env_vars, drop = FALSE]

  # Population SD (matching glmnet / RescaleRasters convention)
  sdN <- function(x) sqrt(mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))

  var_moments <- data.frame(
    Variable = env_vars,
    Mean     = vapply(env_vars, function(v) mean(pm[[v]], na.rm = TRUE), numeric(1)),
    SD       = vapply(env_vars, function(v) sdN(pm[[v]]),                 numeric(1)),
    row.names = NULL
  )

  # Merge moments into coef_tab
  coef_tab <- merge(coef_tab, var_moments, by = "Variable", sort = FALSE)

  # ── 5. Rescale raster layers ─────────────────────────────────────────────────
  pred_subset  <- predictors[[env_vars]]
  scaled_layers <- vector("list", length(env_vars))

  for (i in seq_along(env_vars)) {
    v      <- env_vars[i]
    mu     <- coef_tab$Mean[coef_tab$Variable == v]
    sigma  <- coef_tab$SD[coef_tab$Variable == v]
    beta   <- abs(coef_tab$BetaWeight[coef_tab$Variable == v])

    if (sigma == 0) {
      warning(sprintf("Variable '%s' has zero SD in training data — skipping.", v))
      # Produce a zero layer so downstream clustering ignores it
      scaled_layers[[i]] <- pred_subset[[v]] * 0
    } else {
      # Centre and scale, then weight by |beta|
      scaled_layers[[i]] <- ((pred_subset[[v]] - mu) / sigma) * beta
    }
  }

  pred_rescale        <- terra::rast(scaled_layers)
  names(pred_rescale) <- env_vars

  # ── 6. Optional uncertainty layer ────────────────────────────────────────────
  uncertainty_rast <- NULL
  if (include_uncertainty) {
    uncertainty_rast <- build_uncertainty_layer(
      model        = model,
      pred_mat     = pred_mat,
      predictors   = predictors,
      env_vars     = env_vars,
      pred_rescale = pred_rescale,
      coef_tab     = coef_tab,
      uncertainty_wt = uncertainty_wt
    )
    pred_rescale <- c(pred_rescale, uncertainty_rast)
  }

  # ── 7. Return ─────────────────────────────────────────────────────────────────
  list(
    RescaledPredictors = pred_rescale,
    BetaCoefficients   = coef_tab[, c("Variable", "Estimate", "Est.Error",
                                       "Q2.5", "Q97.5", "BetaWeight")],
    UncertaintyLayer   = uncertainty_rast
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

#' Extract and tidy posterior fixed-effect summaries from a brmsfit
#'
#' Drops the Intercept and all GP / sdgp / lscale parameters so only
#' environmental predictor betas remain.
#'
#' @param model brmsfit
#' @param beta_summary Which summary column to use as BetaWeight
#' @return Data frame: Variable, Estimate, Est.Error, Q2.5, Q97.5, BetaWeight
#' @keywords internal
#' @noRd
extract_posterior_betas <- function(model, beta_summary) {
  # brms::fixef returns matrix: rows = params, cols = Estimate/Est.Error/Q2.5/Q97.5
  fe <- as.data.frame(brms::fixef(model, summary = TRUE))
  fe$Variable <- rownames(fe)
  rownames(fe) <- NULL

  # Drop intercept and any GP-related parameters
  # GP smooth parameters appear as "sgp(...)" or contain "gp(" in the name
  gp_pattern  <- "^(Intercept|sgp|sdgp|lscale)"
  fe <- fe[!grepl(gp_pattern, fe$Variable), , drop = FALSE]

  # brms uses b_ prefix internally; fixef() strips it — but double-check
  # and strip any residual "b_" prefix so names match raster layer names
  fe$Variable <- sub("^b_", "", fe$Variable)

  if (nrow(fe) == 0) {
    stop(
      "No environmental fixed effects found in model. ",
      "Ensure the brmsfit contains named predictor terms."
    )
  }

  fe$BetaWeight <- fe[[beta_summary]]
  fe
}


#' Propagate posterior coefficient uncertainty into a raster layer
#'
#' For each environmental predictor, the posterior SD of its coefficient
#' (Est.Error from fixef) represents how uncertain the model is about that
#' variable's importance. We weight each standardised raster by this SD and
#' combine them in quadrature (sqrt of sum of squares) to produce a single
#' "coefficient uncertainty" surface. This is not a full posterior predictive
#' SD (which would require integrating the GP term too), but it captures
#' environmental-covariate uncertainty efficiently without additional MCMC
#' draws.
#'
#' The layer is then rescaled to match the range of the maximum environmental
#' layer, multiplied by `uncertainty_wt`.
#'
#' @param model brmsfit (unused directly; moments come from coef_tab)
#' @param pred_mat data frame of predictor values
#' @param predictors SpatRaster
#' @param env_vars character vector of env predictor names
#' @param pred_rescale already-scaled SpatRaster (used for range reference)
#' @param coef_tab data frame with Variable, Mean, SD, Est.Error, BetaWeight
#' @param uncertainty_wt numeric weight
#' @return SpatRaster named "coef_uncertainty"
#' @keywords internal
#' @noRd
build_uncertainty_layer <- function(
    model, pred_mat, predictors, env_vars,
    pred_rescale, coef_tab, uncertainty_wt
) {
  # For each variable: (standardised_raster * Est.Error)^2 — then sqrt of sum
  sq_layers <- vector("list", length(env_vars))

  for (i in seq_along(env_vars)) {
    v         <- env_vars[i]
    mu        <- coef_tab$Mean[coef_tab$Variable == v]
    sigma     <- coef_tab$SD[coef_tab$Variable == v]
    beta_sd   <- coef_tab$Est.Error[coef_tab$Variable == v]

    if (sigma == 0) {
      sq_layers[[i]] <- predictors[[v]] * 0
    } else {
      std_rast       <- (predictors[[v]] - mu) / sigma
      sq_layers[[i]] <- (std_rast * beta_sd)^2
    }
  }

  uncertainty_rast <- sqrt(terra::app(terra::rast(sq_layers), sum, na.rm = TRUE))

  # Rescale to match the maximum range across environmental layers, * wt
  env_ranges       <- terra::global(pred_rescale, fun = "range", na.rm = TRUE)
  target_range     <- max(abs(env_ranges$min - env_ranges$max)) * uncertainty_wt

  u_range <- terra::global(uncertainty_rast, fun = "range", na.rm = TRUE)
  u_span  <- as.numeric(u_range$max - u_range$min)

  if (u_span > 0) {
    uncertainty_rast <- (uncertainty_rast - u_range$min) /
      u_span * target_range
  }

  names(uncertainty_rast) <- "coef_uncertainty"
  uncertainty_rast
}
