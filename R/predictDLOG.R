#' Predict NDVI for any DOY based on the fitted DLOG parameters.
#' 
#' The parameters include: 
#'  - mn: winter NDVI
#'  - mx: maximum NDVI
#'  - sos: start-of-season
#'  - rsp: initial slope
#'  - eos: end-of-season
#'  - rau: ending slope
#' 
#' @param params Parameters for prediction function (as listed above).
#' @param t Days since planting. 
#' @return Planting date predictions from input parameters, time after DOS.
#' @export
predictDLOG <- function(params, t) {
  mn = as.numeric(params[1])
  mx =  as.numeric(params[2])
  sos =  as.numeric(params[3])
  rsp =  as.numeric(params[4])
  eos =  as.numeric(params[5])
  rau =  as.numeric(params[6])
  preds = mn + (mx - mn) * (1 / (1 + exp(-rsp * (t - sos))) + 1 / (1 + exp(rau * (t - eos))) - 1)
  return(preds)
}