#' Iterates through all training plots, fits a double logistic, and records the
#' fitted parameters.
#'
#' If the model is not able to be fit to the site data, the site is skipped.
#'
#' @param dlog.data An empty data frame for per-site predictions.
#' @return A data frame of the fitted double logistic function for each site.
#' @export
fit_double_logistic <- function(dlog.data) {
  # Iterate through all training plots.
  for (i in 1:length(trainIDs)) {
    print(i)
    this = train %>% filter(ID == trainIDs[i]) %>% distinct(DOY, .keep_all = TRUE)
    if (nrow(this))
      
      # Fit the model based on Beck's et al. equation.
      mod = try(FitDoubleLog(
        x = this$NDVI,
        t = this$DOY,
        weighting = TRUE,
        plot = FALSE
      ),
      silent = T)
    
    # If we aren't able to fit the model, we try 1x more, then skip if it's
    # still bad.
    if (is.na(mod[2])) {
      mod = try(FitDoubleLog(
        x = this$NDVI,
        t = this$DOY,
        weighting = TRUE,
        plot = FALSE
      ),
      silent = T)
      if (is.na(mod[2])) {
        next
      }
    }
    
    # Calculate the squared error of the fitted values.
    sqerror = (this$NDVI - mod$predicted) ^ 2
    
    # Get the fitted model parameters.
    mod = as.vector(mod$params)
    
    # Get the predicted values for the entire year from the fitted dlog model.
    preds = predictDLOG(mod, c(0:365))
    slope = diff(preds)
    
    # Calculate the departure from the winter NDVI
    # (over-time, to help us determine when the crop likely emerged).
    diffs = preds - as.numeric(mod[1])
    
    # Saves the data for this plot in the 'dlog.data' data frame.
    dlog.data[i,] = data.frame(
      ID = this$ID[1],
      mn =  as.numeric(mod[1]),
      mx =  as.numeric(mod[2]),
      sos =  as.numeric(mod[3]),
      rsp =  as.numeric(mod[4]),
      eos =  as.numeric(mod[5]),
      rau =  as.numeric(mod[6]),
      dos = this$dos[1],
      lat = this$Latitude[1],
      rmse = sqrt(mean((sqerror))),
      maxDOY = which.max(preds),
      emerg = which(slope > 0.0001)[1]
    )
  }
  return(dlog.data)
}
