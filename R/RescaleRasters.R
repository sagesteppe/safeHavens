#' Rescale a raster stack to reflect the beta coefficients from a glmnet model
#' 
#' @description These rescaled rasters can then be used for clustering, and predicting the 
#' results of cluster analysis back into space for a final product. 
#' @param model the final output model from glmnet from `elasticSDM`
#' @param predictors the raster stack to use for the process from `elasticSDM`
#' @param training_data the same data that went into the glmnet model, this is used
#' for calculating variance which is required for the scaling process. From `elasticSDM`
#' @param pred_mat the Prediction matrix from `elasticSDM`
#' @return A list with two objects. 1) The rescaled raster stack. 
#' 2) A table of both standardized and unstandardized coefficients from the glmnet model. 
#' @export 
RescaleRasters <- function(model, predictors, training_data, pred_mat){
  
  sdN <- function(x){
    sqrt((1 / length(x)) * sum((x - mean(x))^2))
  }
  
  # extract coefficients (including intercept)
  coef_mat <- as.matrix(stats::coef(model))

  coef_vec <- coef_mat[, 1]
  coef_names <- rownames(coef_mat)

  coef_tab <- data.frame(
    Variable = coef_names,
    Coefficient = coef_vec,
    row.names = NULL
  )

  # drop intercept explicitly
  coef_tab <- coef_tab[coef_tab$Variable != "(Intercept)", , drop = FALSE]

  # response variance
  yvar <- sdN(as.numeric(training_data$occurrence) - 1)

  # ensure predictors exist
  stopifnot(all(coef_tab$Variable %in% colnames(pred_mat)))

  # predictor SDs (name-safe)
  x_sd <- vapply(
    coef_tab$Variable,
    function(v) sdN(pred_mat[, v]),
    numeric(1)
  )

  # compute beta coefficients (one per row)
  coef_tab$BetaCoefficient <- (coef_tab$Coefficient / yvar) * x_sd

  # this rescales the raster to be equivalent to the inputs to the elastic net model. 
  # after this they still need to be multiplied by the beta coefficients 
  stopifnot(inherits(predictors, "SpatRaster"))

  pred_rescale <- terra::subset(
    predictors,
    coef_tab$Variable
  )

  layer_names <- names(pred_rescale)
  out_layers <- lapply(layer_names, function(lyr_name) {

    vals <- pred_mat[, lyr_name]

    scaled <- terra::app(
      pred_rescale[[lyr_name]],
      fun = function(x){ (x - mean(vals)) / sdN(vals) }
    )

    beta <- abs(
      coef_tab$BetaCoefficient[
        coef_tab$Variable == lyr_name
      ]
    )

    scaled * beta
  })
  pred_rescale <- terra::rast(out_layers)
  names(pred_rescale) <- layer_names


  list(
    RescaledPredictors = pred_rescale, 
    BetaCoefficients = coef_tab
  )
}

