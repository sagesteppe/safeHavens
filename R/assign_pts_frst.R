#' Assign randomly sampled points to their nearest grid based on the difference
#' between the polygons original and desired size. 
#' 
#' @param x distance matrix between the points and their potential neighbors polygons
#' @param props the existing proportional area of each potential neighbor polygons
#' @param nf_pct the desired proportional area of each potential neighbor polygon
#' @keywords internal
#' @noRd
assign_pts_frst <- function(x, props, nf_pct){

  # ensure these are in the same order. so we can match by position
  # in the loop. 
  need_most2least <- names(sort(nf_pct - props))
  props <- props[need_most2least]
  x <- x[,need_most2least]
  
  x$Assignment  <- NA
  x$ID <- seq_len(nrow(x))

  for(i in seq_along(props)){
    
    # we assign the grids to the neediest grids first, and then work back 
    # removing these points so they are not overwritten. 
    x_sub <- x[is.na(x$Assignment),]

        # ensure we don't exceed available points
    n_assign <- min(props[i], nrow(x_sub))
    if(n_assign > 0){
      indices <- x_sub[sort(x_sub[, i], index.return = TRUE)$ix[1:n_assign], 'ID']
      x[indices, 'Assignment'] <- names(props)[i]
    }
  }
  x
}
