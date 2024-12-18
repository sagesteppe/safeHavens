#' Save the results of SDMs from the `elasticSDM`, `RescaleRasters` and `PostProcessSDM` function in safeHavens
#' 
#' @description This function is used to write out a wide range of values from the `fitPredictOperationalize` process.
#' It will create multiple subdirectories within a user specified path. 
#' These include: 'Rasters' where a raster stack of the four final rasters will go, 'Fitting' where the details of model fitting from caret will be placed, 'Models' where the final fit model will go, 'Evaluation' where all evaluation statistics will be placed, 'Threshold' where results form dismo::threshold will be placed. 
#' @param path a root path where each of 5 folders will be created, if they do not exist. 
#' @param taxon the name of the taxonomic entity for which the models were created. 
#' @param cv_model the cross validation data from `elasticSDM`
#' @param pcnm the pcnm/mem rasters from `elasticSDM`
#' @param model the final glmnet model from `elasticSDM`
#' @param cm the confusion matrix from `elasticSDM`
#' @param coef_tab the coefficient table from `RescaleRasters`
#' @param f_rasts the final rasters from `RescaleRasters`
#' @param thresh threshold statistics from `PostProcessSDM`
#' @return all the above objects, or all objects above specified, are written to disk. 
#' @export
writeSDMresults <- function(path, taxon, cv_model, pcnm, model, cm, coef_tab, f_rasts, thresh){
  
  if(!missing(cv_model)){
    dir.create(file.path(path, 'Fitting'), showWarnings = FALSE)
    saveRDS( # all cross validation information
      cv_model, 
      file = file.path(path, 'Fitting', paste0(taxon, '-Fitting.rds'))   
    ) 
  } else {message('param: `cv_model` not supplied to fn, not saving it.')}
  
  if(!missing(pcnm)){
    dir.create(file.path(path, 'PCNM'), showWarnings = FALSE)
    terra::writeRaster( # save the pcnm layers we created. 
      pcnm, overwrite = TRUE, 
      file = file.path(path, 'PCNM', paste0(taxon, '-PCNM.tif'))
    ) 
  } else {message('param: `pcnm` not supplied to fn, not saving it.')}

  if(!missing(model)){
    dir.create(file.path(path, 'Model'), showWarnings = FALSE)
    saveRDS( # save the final fitted model.
      model, 
      file = file.path(path, 'Model', paste0(taxon, '-Model.rds'))
    )
  } else {message('param: `model` not supplied to fn, not saving it.')}
  
  if(!missing(coef_tab)){ # Model coefficients 
    dir.create(file.path(path, 'Model'), showWarnings = FALSE)
    utils::write.csv( # save the model coefficients
      coef_tab,  row.names = FALSE,
      file = file.path(path, 'Model', paste0(taxon, '-Coefficients.csv'))
    )
  } else {message('param: `coef_tab` not supplied to fn, not saving it.')}

  if(!missing(cm)){ # Confusion Matrix Materials here. 
    
    dir.create(file.path(path, 'Evaluation'), showWarnings = FALSE)
    
    utils::write.csv( #  save confusion matrix from old school split. 
      data.frame(t(cm$table)), row.names = FALSE, 
      file = file.path(path, 'Evaluation', paste0(taxon, '-CMatrix.csv'))
    )
      
    utils::write.csv( # save the calculated evaluation metrics. 
      data.frame(
        Variable = names(cm$byClass),
        Value = as.numeric(cm$byClass)
      ), row.names = FALSE, 
      file = file.path(path, 'Evaluation', paste0(taxon, '-Metrics.csv'))
    )  
  }  else {message('param: `cm` not supplied to fn, not saving it.')}
  
  if(!missing(thresh)){
    dir.create(file.path(path, 'Threshold'), showWarnings = FALSE)
    utils::write.csv( # threshold information - make it tidy. 
      data.frame(
        Metric = colnames(thresh), 
        Value = as.numeric(t(thresh))
      ), row.names = FALSE,
      file = file.path(path, 'Threshold', paste0(taxon, '-Threshold.csv'))
    )
  } else {message('param: `thresh` not supplied to fn, not saving it.')}
  
  if(!missing(f_rasts)){
    dir.create(file.path(path, 'Raster'), showWarnings = FALSE)
    terra::writeRaster( # save all final rasters
      f_rasts, overwrite = TRUE, 
      filename = file.path(path, 'Raster', paste0(taxon, '.tif'))
    )
  } else {message('param: `f_rasts` not supplied to fn, not saving it.')}
  
}

