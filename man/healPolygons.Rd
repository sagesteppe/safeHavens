% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/healPolygons.R
\name{healPolygons}
\alias{healPolygons}
\title{Clean up unioned geometries - part 1}
\usage{
healPolygons(x)
}
\arguments{
\item{x}{most of the output of \code{snapgrids}}
}
\description{
this function uses sf::st_snap to remove small lines and other artifacts associated
with the unioning of polygons. This is ran within \code{snapGrids}
}
\keyword{internal}
