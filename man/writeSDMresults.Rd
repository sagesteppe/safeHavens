% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WriteSDMresults.R
\name{writeSDMresults}
\alias{writeSDMresults}
\title{Save the results of SDMs from the \code{elasticSDM}, \code{RescaleRasters} and \code{PostProcessSDM} function in safeHavens}
\usage{
writeSDMresults(
  path,
  taxon,
  cv_model,
  pcnm,
  model,
  cm,
  coef_tab,
  f_rasts,
  thresh
)
}
\arguments{
\item{path}{a root path where each of 5 folders will be created, if they do not exist.}

\item{taxon}{the name of the taxonomic entity for which the models were created.}

\item{cv_model}{the cross validation data from \code{elasticSDM}}

\item{pcnm}{the pcnm/mem rasters from \code{elasticSDM}}

\item{model}{the final glmnet model from \code{elasticSDM}}

\item{cm}{the confusion matrix from \code{elasticSDM}}

\item{coef_tab}{the coefficient table from \code{RescaleRasters}}

\item{f_rasts}{the final rasters from \code{RescaleRasters}}

\item{thresh}{threshold statistics from \code{PostProcessSDM}}
}
\value{
all the above objects, or all objects above specified, are written to disk.
}
\description{
This function is used to write out a wide range of values from the \code{fitPredictOperationalize} process.
It will create multiple subdirectories within a user specified path.
These include: 'Rasters' where a raster stack of the four final rasters will go, 'Fitting' where the details of model fitting from caret will be placed, 'Models' where the final fit model will go, 'Evaluation' where all evaluation statistics will be placed, 'Threshold' where results form dismo::threshold will be placed.
}
