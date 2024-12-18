% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnvironmentalBasedSample.R
\name{EnvironmentalBasedSample}
\alias{EnvironmentalBasedSample}
\title{Create environmental and spatial clusters for targeting collection areas}
\usage{
EnvironmentalBasedSample(
  pred_rescale,
  f_rasts,
  lyr,
  taxon,
  path,
  n,
  fixedClusters,
  n_pts,
  planar_projection,
  coord_wt,
  buffer_d,
  prop_split,
  write2disk,
  ...
)
}
\arguments{
\item{pred_rescale}{a rasterstack of predictor layers which have been rescaled to represent the beta coefficients from the elastic net (glmnet::glmnet) modelling process. See ?RescaleRasters for an implementation of this functionality.}

\item{f_rasts}{The rasters output from the SDM workflow.}

\item{lyr}{Character. The name of the layer you want to use for the analysis one of: 'Threshold', 'Clipped', or Supplemented'. If missing defaults to 'Supplemented'.}

\item{taxon}{Character. the name of the taxonomic entity for which the models were created. The final raster of clusters, the results from both KNN classifier trainings, and details of the clustering procedure (if fixedClusters=TRUE).}

\item{path}{a root path where each output data will be saved, use the same as in \code{WriteSDMresults}. Defaults to current working directory, check this with getwd().}

\item{n}{Numeric. the number of clusters desired.}

\item{fixedClusters}{Boolean. Defaults to TRUE, which will create n clusters. If False then use NbClust::NbClust to determine the number of clusters.}

\item{n_pts}{Numeric. the number of points to use for generating the clusters, these will be randomly sampled within the mask area \code{mask}. Defaults to 500.}

\item{planar_projection}{Numeric, or character vector. An EPSG code, or a proj4 string, for a planar coordinate projection, in meters, for use with the function. For species with very narrow ranges a UTM zone may be best (e.g. 32611 for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). Otherwise a continental scale projection like 5070 See https://projectionwizard.org/ for more information on CRS. The value is simply passed to sf::st_transform if you need to experiment.}

\item{coord_wt}{Numeric. The amount to weigh coordinates by for the distance matrix relative to the most important variable identified by elastic net regression. Defaults to 2.5, setting at 1 would make this value equivalent to environmental and PCNM variables. This metric increases the spatial contiguity of the clusters identified.}

\item{buffer_d}{Numeric. When using two-stage sampling to increase the sample size (number of points) for uncommon clusters, the distance to buffer the individual pts which are located in these cells for further possible sampling. Defaults to 3 which allows for each an area encompassing around 45-50 raster cells nearby to possibly be sampled. This is a reasonable default for coarsely gridded data (the most sensible grain type for this packages use cases), moderate resolution data (e.g. 250m-1km) may require higher values to find meaningful differences between the variables at these further locations.}

\item{prop_split}{Numeric. The proportion of records to be used for training the KNN classifier. Defaults to
0.8 to use 80\% of records for training and 20\% for the independent test sample.}

\item{write2disk}{Boolean. Whether to write results to disk or not. Defaults to FALSE.}

\item{...}{Further arguments passed to NbClust::NbClust for optimizing cluster numbers. Defaults to using method = 'complete', which will compare the results from 20 methods and select the cluster number most commonly generated by all of these algorithms. We have min.nc set as a default of 5, to overcome having too few clusters to use (but this easily overwritten by supplying the argument minc.nc=2, to set it lower), and max.nc = 20 for congruence with many seed collection endeavors (overwritten by max.nc = 10 for example).}
}
\value{
Writes four objects to disk, and returns one object to R session (optional).
}
\description{
This function utilizes the output from an elastic net GLM model to create a weights matrix of
features relevant to a species distribution to identify clusters throughout it's range while
incorporating the PCNM/MEM data and coordinates to implement some spatial contiguity.
}
