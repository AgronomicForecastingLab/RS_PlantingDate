#' Applies a uniform 3x3 smoothing window across the input matrix
#' (excluding pixels equal to NA or -999).
#'
#' A smaller smoothing window (e.g., 2x2, 2x3) is applied to edge and corner pixels.
#'
#' @param dos.mat The matrix of pixels to be smoothed. 
#' @return A smoothed matrix (excluding NA/-999 pixels).
#' @export
smooth_3x3 <- function(dos.mat) {
  smooth.mat = dos.mat * 0
  na.mat = dos.mat
  na.mat[which(dos.mat == -999)] = NA
  
  # Store dimensions of the input matrix.
  d1 = nrow(dos.mat)
  d2 = ncol(dos.mat)
  
  # Smooths the matrix with a 3x3 window, skipping corners & edges.
  for (r in 2:(d1 - 1)) {
    print(r)
    for (c in 2:(d2 - 1)) {
      # Skips smoothing if the pixel is NA or -999.
      if (is.na(dos.mat[r, c]))
        next
      if (dos.mat[r, c] == -999)
        next
      
      # Take the average of the surrounding 3x3 window.
      temp = na.mat[(r - 1):(r + 1), (c - 1):(c + 1)]
      smooth.mat[r, c] = round(mean(temp, na.rm = TRUE))
    }
  }
  
  # Smoothes the edges of the prediction matrix (reducing the 3x3 window to an
  # appropriately sized window).
  # Left edges:
  for (r in 2:(d1 - 1)) {
    if (is.na(dos.mat[r, 1]))
      next
    if (dos.mat[r, 1] == -999)
      next
    
    # Get 2x3 surrounding matrix and take the average.
    temp = na.mat[(r - 1):(r + 1), 1:2]
    smooth.mat[r, 1] = round(mean(temp, na.rm = TRUE))
  }
  
  # Right edges:
  for (r in 2:(d1 - 1)) {
    if (is.na(dos.mat[r, d2]))
      next
    if (dos.mat[r, d2] == -999)
      next
    
    # Get 2x3 surrounding matrix and take the average.
    temp = na.mat[(r - 1):(r + 1), (d2 - 1):d2]
    smooth.mat[r, d2] = round(mean(temp, na.rm = TRUE))
  }
  
  # Top edges:
  for (c in 2:(d2 - 1)) {
    if (is.na(dos.mat[1, c]))
      next
    if (dos.mat[1, c] == -999)
      next
    
    # Get 3x2 surrounding matrix and take the average.
    temp = na.mat[1:2, (c - 1):(c + 1)]
    smooth.mat[1, c] = round(mean(temp, na.rm = TRUE))
  }
  
  # Bottom edges:
  for (c in 2:(d2 - 1)) {
    if (is.na(dos.mat[d1, c]))
      next
    if (dos.mat[d1, c] == -999)
      next
    
    # Get 3x2 surrounding matrix and take the average.
    temp = na.mat[(d1 - 1):d1, (c - 1):(c + 1)]
    smooth.mat[d1, c] = round(mean(temp, na.rm = TRUE))
  }
  
  # Smoothes the corners of the prediction matrix (reducing the 3x3 window to a
  # 2x2 matrix).
  # Top left corner:
  if (!is.na(dos.mat[1, 1]) & dos.mat[1, 1] != -999) {
    temp = na.mat[1:(r + 1), 1:(c + 1)]
    smooth.mat[1, 1] = round(mean(temp, na.rm = TRUE))
  }
  # Top right corner:
  if (!is.na(dos.mat[1, d2]) & dos.mat[1, d2] != -999) {
    temp = na.mat[1:(r + 1), (d2 - 1):d2]
    smooth.mat[1, d2] = round(mean(temp, na.rm = TRUE))
  }
  # Bottom left corner:
  if (!is.na(dos.mat[d1, 1]) & dos.mat[d1, 1] != -999) {
    temp = na.mat[(d1 - 1):d1, 1:(c + 1)]
    smooth.mat[d1, 1] = round(mean(temp, na.rm = TRUE))
  }
  # Bottom right corner:
  if (!is.na(dos.mat[d1, d2]) & dos.mat[d1, d2] != -999) {
    temp = na.mat[(d1 - 1):d1, (d2 - 1):d2]
    smooth.mat[d1, d2] = round(mean(temp, na.rm = TRUE))
  }
  return(smooth.mat)
}
