#' Create a simple, theoretical, raster surface modelling Isolation by Distance.
#'
#' Used in conjunction, preferably before `populationResistance` in the `IBRBasedSample` workflow.
#'
#' @param base_raster SpatRaster. Base raster for the study area. Provides template geometry and resolution.
#' @param resistance_surface SpatRaster. Optional pre-computed resistance raster. If provided, the raster-building arguments are ignored.
#' @param oceans SpatRaster. Binary (0/1) raster for ocean cells. Used to increase movement cost.
#' @param lakes SpatRaster. Binary (0/1) raster for lakes.
#' @param rivers SpatRaster. Binary (0/1) raster for rivers.
#' @param tri SpatRaster. Continuous raster of topographic roughness (TRI). Used to increase cost in mountainous terrain.
#' @param habitat SpatRaster. Continuous raster of habitat suitability. Low values increase cost.
#' @param addtl_r SpatRaster, 'raster stack'. Additional layers to include in the resistance surface
#' @param w_ocean Numeric. Weight applied to oceans (default 2000).
#' @param w_lakes Numeric. Weight applied to lakes (default 200).
#' @param w_rivers Numeric. Weight applied to rivers (default 20).
#' @param w_tri Numeric. Weight applied to TRI (default 1).
#' @param w_habitat Numeric. Weight applied to habitat suitability (default 1).
#' @param addtl_w Numeric vector. Must equal the length of addtl_r exactly. Weights for the additional rasters layers to include in the resistance surface

#' @examples
#' \dontrun{
#' # Prepare resistance raster
#'
#' # this also can run internally in `population resistance`,
#' # but for time sakes is best to prep ahead of time
#' # especially if treating multiple species in the same domain.
#' res <- buildResistanceSurface(
#'   base_raster = base_rast,
#'   oceans = ocean_r,
#'   lakes = lakes_r,
#'   rivers = rivers_r,
#'   tri = tri_r
#' )
#' }
#' @export
buildResistanceSurface <- function(
  base_raster,
  resistance_surface = NULL,
  oceans = NULL,
  lakes = NULL,
  rivers = NULL,
  tri = NULL,
  habitat = NULL,
  addtl_r = NULL,
  w_ocean = 100,
  w_lakes = 50,
  w_rivers = 20,
  w_tri = 1,
  w_habitat = 1,
  addtl_w = NULL
) {
  if (!is.null(resistance_surface)) {
    res <- resistance_surface

    if (!terra::compareGeom(res, base_raster, stopOnError = FALSE)) {
      stop("resistance_surface must match base_raster geometry")
    }
  } else {
    # initialize blank raster
    res <- terra::rast(base_raster)
    terra::values(res) <- 0

    # add weighted features
    if (!is.null(oceans)) {
      res[oceans == 1] <- res[oceans == 1] + w_ocean
    }
    if (!is.null(lakes)) {
      res[lakes == 1] <- res[lakes == 1] + w_lakes
    }
    if (!is.null(rivers)) {
      res[rivers == 1] <- res[rivers == 1] + w_rivers
    }
    if (!is.null(tri)) {
      res <- res + w_tri * terra::setValues(res, scale(terra::values(tri)))
    }
    if (!is.null(habitat)) res <- res + w_habitat * habitat
    
    ## additional layers if supplied
    if (!is.null(addtl_r) & terra::nlyr(addtl_r)==length(addtl_w)){
      for (i in seq_along(addtl_w)){
        res <- res + addtl_r[[i]] * addtl_w[i]
      }
    }
  }

  # --- clamp only if raster has values ---
  vals <- terra::values(res)
  if (!is.null(vals) && length(vals) > 0 && any(!is.na(vals))) {
    res <- terra::clamp(res, lower = 1L)
  }

  # convert to integer for memory and calculation efficiency.
  res <- terra::as.int(res)
}
