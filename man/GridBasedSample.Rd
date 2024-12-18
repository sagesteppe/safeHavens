% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GridBasedSample.R
\name{GridBasedSample}
\alias{GridBasedSample}
\title{Create hexagonal grid based polygons over a geographic range}
\usage{
GridBasedSample(x, planar_projection, gridDimensions)
}
\arguments{
\item{x}{An SF object or terra spatraster. the range over which to generate the clusters.}

\item{planar_projection}{Numeric, or character vector. An EPSG code, or a proj4 string, for a planar coordinate projection, in meters, for use with the function. For species with very narrow ranges a UTM zone may be best (e.g. 32611 for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). Otherwise a continental scale projection like 5070 See https://projectionwizard.org/ for more information on CRS. The value is simply passed to sf::st_transform if you need to experiment.}

\item{gridDimensions}{A single row form the ouput of \code{TestGridSizes} with the optimal number of grids to generate.}
}
\value{
An simple features (sf) object containing the final grids for saving to computer. See the vignette for questions about saving the two main types of spatial data models (vector - used here, and raster).
}
\description{
This function creates 20 grid cells over a geographic area (\code{x}), typically a species range.
}
\examples{
ri <- spData::us_states |> 
dplyr::filter(NAME == 'Rhode Island') |>
  sf::st_transform(32615)

sizeOptions <- TestGridSizes(ri)
head(sizeOptions) # in this case let's shoot for 33 and see what happens
sizeOptions <- sizeOptions[sizeOptions$Name == 'Original',]

output <- GridBasedSample(ri, 5070, gridDimensions = sizeOptions)
plot(output)
}
