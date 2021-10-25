#' Parse through all pixels and record optimized DOS.
#' 
#' Calculates optimized DOS by averaging the results of both Venancio
#' equations. Generates predictions for pixels with at least 3 dates of data. 
#' 
#' Predictions outside of a reasonable DOS range (before/after insurance 
#' coverage dates) are set to -999. 
#'
#' @param dat A x-y-n array of NDVI data. 
#' @param doys A list of days since planting at each observed date.
#' @param fA Function for 1st Venancio equation.
#' @param fB Function for 2nd Venancio equation.
#' @param start_d The earliest DOS value for the observed region.
#' @param end_d The latest DOS value for the observed region.
#' @return A matrix of DOS predictions. 
#' @export
predict_DOS  <- function(dat, doys, fA, fB, start_d, end_d) {
  dos.mat = matrix(NA, dim(dat)[1], dim(dat)[2])
  for (i in 1:dim(dat)[1]) {
    print(i)
    for (j in 1:dim(dat)[2]) {
      # Get data for this pixel.
      pixel = data.frame(t = as.numeric(doys), ndvi = dat[i, j, ])
      
      # Skip over masked pixels (pixels with fewer than 3 valid data points)
      if (length(which(is.na(pixel$ndvi))) >= 21) {
        next
      }
      
      # Set the bare-ground value to the first observed date for the pixel.
      if (!is.na(pixel$ndvi[1])) {
        preNDVI = pixel$ndvi[1]
      } else {
        preNDVI = pixel$ndvi[2]
      }
      
      # Optimize the RMSE for the two functions and record the optimized value &
      # associated RMSE.
      opA = optimize(fA,
                     c(0, 260),
                     t = pixel$t,
                     N = pixel$ndvi,
                     minNDVI = preNDVI)
      opB = optimize(fB,
                     c(0, 260),
                     t = pixel$t,
                     N = pixel$ndvi,
                     minNDVI = preNDVI)
      
      # Corn typically isn't planted before 04/01 or after 06/06 in IL.
      # Set any early values to NA. Otherwise, set value.
      outside = (round(mean(opA$minimum, opB$minimum)) < start_e) |
        (round(mean(opA$minimum, opB$minimum)) > end_d)
      dos.mat[i, j] = ifelse(outside,-999, round(mean(opA$minimum, opB$minimum)))
    }
  }
  return(dos.mat)
}
