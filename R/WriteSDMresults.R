#' Save the results of SDMs from the 'fitPredictOperationalize' function in safeHavens
#' 
#' This function is used to write out a wide range of values from the 
#' `fitPredictOperationalize` process. It will create multiple subdirectories 
#' within a user specified path. 
#' These include: 'Rasters' where a raster stack of the four final rasters will go, 
#' Fitting' where the details of model fitting from Caret will be placed, 
#' 'Models' where the final fit model will go, 'Evaluation' where all evaluation 
#' statistics will be placed, 'Threshold' where results form dismo::threshold will
#' be placed. 
#' @param path a root path where each of 5 folders will be created, if they do not exist. 
#' @param taxon the name of the taxonomic entity for which the models were created. 
#' 
writeSDMresults <- function(path, taxon){
  
  dir.create(file.path(path, 'PCNM'), showWarnings = FALSE)
  dir.create(file.path(path, 'Fitting'), showWarnings = FALSE)
  dir.create(file.path(path, 'Model'), showWarnings = FALSE)
  dir.create(file.path(path, 'Evaluation'), showWarnings = FALSE)
  dir.create(file.path(path, 'Threshold'), showWarnings = FALSE)
  dir.create(file.path(path, 'Raster'), showWarnings = FALSE)
  
  saveRDS( # all cross validation information
    cv_model, 
    file = file.path(path, 'Fitting', paste0(taxon, '-Fitting.rds'))   
    ) 
  
  terra::writeRaster( # save the pcnm layers we created. 
    pcnm, overwrite = TRUE, 
    file = file.path(path, 'PCNM', paste0(taxon, '-PCNM.tif'))
    ) 
  
  saveRDS( # save the final fitted model.
    mod, 
    file = file.path(path, 'Model', paste0(taxon, '-Model.rds'))
    )
  
  write.csv( # save the model coefficients
    coef_tab,  row.names = FALSE,
    file = file.path(path, 'Model', paste0(taxon, '-Coefficients.csv'))
  )
  
  write.csv( #  save confusion matrix from old school split. 
    data.frame(t(cm$table)), row.names = FALSE, 
    file = file.path(path, 'Evaluation', paste0(taxon, '-CMatrix.csv'))
  ) 
  
  write.csv( # save the calculated evaluation metrics. 
    data.frame(
      Variable = names(cm$byClass),
      Value = as.numeric(cm$byClass)
    ), row.names = FALSE, 
    file = file.path(path, 'Evaluation', paste0(taxon, '-Metrics.csv'))
  )
  
  write.csv( # threshold information - make it tidy. 
    data.frame(
      Metric = colnames(thresh), 
      Value = as.numeric(t(thresh))
      ), row.names = FALSE,
    file = file.path(path, 'Threshold', paste0(taxon, '-Threshold.csv'))
  )
  
  terra::writeRaster( # save all final rasters
    f_rasts, overwrite = TRUE, 
    filename = file.path(path, 'Raster', paste0(taxon, '.tif'))
  )
}

