#' Iterates through all training plots, fits a normal dist, and records the
#' fitted parameters.
#'
#' If the model is not able to be fit to the site data, the site is skipped.
#'
#' @param x vector of NDVI values
#' @param t vector of DOY values with NDVI values
#' @return A data frame of the fitted normal function for each site.
#' @export
fit_normal <- function(x, t) {
  
  # take a guess at what mu will be (peak NDVI DOY); we use 200 if it's too low 
  #estMu = ifelse(this$DOY[which.max(this$NDVI)] > 172,this$DOY[which.max(this$NDVI)], 200)
  estMu = t[which.max(x)]
  dat = data.frame(NDVI = x, DOY = t)
  res <- try(nls(NDVI ~ k*exp(-1/2*(DOY-mu)^2/sigma^2), 
                 start=c(mu=estMu,sigma=30,k=1.5),
                 nls.control(maxiter = 1500), 
                 data = dat), 
             silent = TRUE)
  if (is.na(res[2])){
    return(NA)
  }
  
  # get fitted model parameters
  res = data.frame(ID = this$ID[1],
                   mu = summary(res)$parameters[1,1],
                   sigma = summary(res)$parameters[2,1],
                   k = summary(res)$parameters[3,1],
                   dos = this$dos[1],
                   lat = this$Latitude[1])
  
  return(res)
}


