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
  
  sdN <- function(x){sigma=sqrt((1/length(x)) * sum((x-mean(x))^2))}
  
  coef_tab <- data.frame(
    Variable = row.names(as.data.frame(as.matrix(stats::coef(model)))), 
    Coefficient = as.numeric(stats::coef(model))
  )
  coef_tab <- coef_tab[2:nrow(coef_tab),]
  
  # now calculate the beta coefficient 
  yvar <- sdN(as.numeric(training_data$occurrence)-1)
  coef_tab$BetaCoefficient <- apply(pred_mat, 2, FUN = sdN)
  coef_tab$BetaCoefficient <- coef_tab$Coefficient / yvar * coef_tab$BetaCoefficient
  
  # this rescales the raster to be equivalent to the inputs to the elastic net model. 
  # after this they still need to be multiplied by the beta coefficients 
  pred_rescale <- predictors
  pred_rescale <- pred_rescale[[ names(pred_rescale) %in% coef_tab$Variable ]]
  for (i in seq_along(1:dim(pred_rescale)[3])){
    
    lyr_name <- names(pred_rescale)[[i]]
    vals <- pred_mat[,lyr_name]
    pred_rescale[[i]] <- terra::app(pred_rescale[[i]], 
                                    fun = function(x){(x - mean(vals)) / sdN(vals)})
    
    pred_rescale[[i]] <- pred_rescale[[i]] * abs(coef_tab[coef_tab$Variable==lyr_name,'BetaCoefficient'])
    names(pred_rescale[[i]]) <- lyr_name
  }
  return(list(
    RescaledPredictors = pred_rescale, 
    BetaCoefficients = coef_tab)
    )
  
}
