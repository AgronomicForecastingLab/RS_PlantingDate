#' Applies a uniform 3x3 smoothing window across the input matrix
#' (excluding pixels equal to NA or -999).
#'
#' A smaller smoothing window (e.g., 2x2, 2x3) is applied to edge and corner pixels.
#'
#' @param x x-dimension (width) of the NDVI images.
#' @param y y-dimension(length) of the NDVI images.
#' @param n number of observed dates.
#' @return A x-y-n matrix of NDVI data.
#' @export
aggregate_dates <- function(x, y, n) {
  # Read in all the adequate NDVI images into one array.
  dat = array(-999, dim = c(x, y, n))
  
  # Aggregate data from all 24 dates into a matrix
  for (i in 1:n) {
    print(i)
    f = stack(files[i])
    mat = as.matrix(f)
    dat[, , i] = matrix(
      data = mat[, 1],
      byrow = T,
      nrow = x,
      ncol = y,
    )
  }
  return(dat)
}