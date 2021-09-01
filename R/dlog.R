#' Calculates RMSE values for the optimized double logistic function.
#'
#' I.e., evaluates the success of the double-logistic function at predicting
#' DOS.
#'
#' #' The parameters include:
#'  1) mn: winter NDVI
#'  2) mx: maximum NDVI
#'  3) sos: start-of-season
#'  4) rsp: initial slope
#'  5) eos: end-of-season
#'  6) rau: ending slope
#'
#' @param par Parameters for prediction function (as listed above).
#' @param t Days since planting.
#' @param NDVI Actual observed NDVI for the site at a time t.
#' @return The RMSE of the double-logistic function predictions.
#' @export
dlog <- function(par, t, NDVI) {
  preds = par[1] + (par[2] - par[1]) * (1 / (1 + exp(-par[4] * (t - par[3]))) + 1 /
                                          (1 + exp(par[6] * (t - par[5]))) - 1)
  sqerror = (preds - NDVI) ^ 2
  return(sqrt(mean(sqerror)))
}