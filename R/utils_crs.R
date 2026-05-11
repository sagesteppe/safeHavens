# Normalise a user-supplied CRS value for sf or terra.
#
# sf accepts bare numeric EPSG codes and prefixed strings ("EPSG:", "ESRI:").
# terra requires a prefixed string; bare integers must become "epsg:<n>".

#' @keywords internal
#' @noRd
planar_proj_sf <- function(x) x

#' @keywords internal
#' @noRd
planar_proj_terra <- function(x) {
  if (is.numeric(x)) paste0("epsg:", x) else x
}
