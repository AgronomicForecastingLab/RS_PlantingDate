require(plyr)

#' Fits a double logistic curve to observed values using a function described in
#' Beck et al. (2006) (equation 3).
#'
#' @param x The vector (or time-series) to fit.
#' @param t Time steps.
#' @param tout Time steps of output (can be used for interpolation).
#' @param weight Weighting scheme (Beck), higher values get higher weight, etc.
#' @param hessian Hessian matrix used to compute standard error of parameters.
#' @param plot If true, the plot is iterating until best logistic fit.
#' @param ninit Number of initial params on which to optimize DLOG.
#' @start_d The earliest DOS value for the observed region.
#' @end_d The latest DOS value for the observed region.
#' @return A matrix of DOS predictions.
#' @export
optimize_double_log <-
  function(x,
           t,
           tout = t,
           weight = TRUE,
           hessian = FALSE,
           plot = FALSE,
           ninit = 30) {
    n <- length(x)
    avg <- mean(x, na.rm = TRUE)
    mx <- max(x, na.rm = TRUE)
    mn <- min(x, na.rm = TRUE)
    ampl <- mx - mn
    .doubleLog <- function(par, t) {
      mn <- par[1]
      mx <- par[2]
      sos <- par[3]
      rsp <- par[4]
      eos <- par[5]
      rau <- par[6]
      xpred <- mn + (mx - mn) * (1 / (1 + exp(-rsp * (t - sos))) +
                                   1 / (1 + exp(rau * (t - eos))) - 1)
      return(xpred)
    }
    .error <- function(par, x, weights) {
      if (any(is.infinite(par)))
        return(99999)
      if (par[1] > par[2])
        return(99999)
      xpred <- .doubleLog(par, t = t)
      sse <- sum((xpred - x) ^ 2 * weights, na.rm = TRUE)
      return(sse)
    }
    if (weight) {
      iter <- 1:2
    } else {
      iter <- 1
    }
    weights <- rep(1, length(x))
    doy <- quantile(t, c(0.25, 0.75), na.rm = TRUE)
    prior <- rbind(
      c(mn, mx, doy[1], 0.5, doy[2], 0.5),
      c(mn,
        mx, doy[2], 0.5, doy[1], 0.5),
      c(mn - ampl / 2, mx + ampl / 2,
        doy[1], 0.5, doy[2], 0.5),
      c(mn - ampl / 2, mx + ampl / 2,
        doy[2], 0.5, doy[1], 0.5)
    )
    prior <- rbind(prior,
                   apply(prior, 2, function(x)
                     rnorm(ninit - 4, mean(x), sd(x) * 2)))
    if (plot)
      plot(t, x)
    for (i in iter) {
      if (i > 1) {
        opt.df <- opt.df[order(opt.df$cost), ]
        prior <- rbind(opt.df[1:5, -(1:2)], opt$par)
      }
      opt.l <-
        apply(
          prior,
          1,
          optim,
          .error,
          x = x,
          weights = weights,
          method = "L-BFGS-B",
          lower = c(0, 0, 0, 0, 0, 0),
          upper = c(1, 1, 365, Inf, 365, Inf),
          control = list(maxit = 2500),
          hessian = FALSE
        )
      opt.df <-
        cbind(
          cost = unlist(llply(opt.l, function(opt)
            opt$value)),
          convergence = unlist(llply(opt.l, function(opt)
            opt$convergence)),
          ldply(opt.l, function(opt)
            opt$par)
        )
      best <- which.min(opt.df$cost)
      if (opt.df$convergence[best] == 1) {
        opt <- opt.l[[best]]
        opt <-
          optim(
            opt.l[[best]]$par,
            .error,
            x = x,
            weights = weights,
            method = "L-BFGS-B",
            control = list(maxit = 1500),
            hessian = FALSE
          )
        prior <- rbind(prior, opt$par)
        xpred <- .doubleLog(opt$par, t)
      } else if (opt.df$convergence[best] == 0) {
        opt <- opt.l[[best]]
        prior <- rbind(prior, opt$par)
        xpred <- .doubleLog(opt$par, t)
      }
      parinit <- opt$par
      mn <- opt$par[1]
      mx <- opt$par[2]
      sos <- opt$par[3]
      rsp <- opt$par[4]
      eos <- opt$par[5]
      rau <- opt$par[6]
      m <- lm(c(0, 100) ~ c(sos, eos))
      tr <- coef(m)[2] * t + coef(m)[1]
      tr[tr < 0] <- 0
      tr[tr > 100] <- 100
      res <- xpred - x
      weights <- 1 / ((tr * res + 1) ^ 2)
      weights[res > 0 & res <= 0.01] <- 1
      weights[res < 0] <- 4
    }
    if (hessian) {
      opt <- optim(
        opt$par,
        .error,
        x = x,
        weights = weights,
        method = "BFGS",
        hessian = TRUE
      )
      .qr.solve <- function(a,
                            b,
                            tol = 1e-07,
                            LAPACK = TRUE) {
        if (!is.qr(a))
          a <- qr(a, tol = tol, LAPACK = LAPACK)
        nc <- ncol(a$qr)
        nr <- nrow(a$qr)
        if (a$rank != min(nc, nr))
          stop("singular matrix 'a' in solve")
        if (missing(b)) {
          if (nc != nr)
            stop("only square matrices can be inverted")
          b <- diag(1, nc)
        }
        res <- try(qr.coef(a, b), silent = TRUE)
        if (class(res)[1] == "try-error") {
          res <- matrix(NA, ncol(b), nrow(b))
        }
        res
      }
      vc <- .qr.solve(opt$hessian)
      npar <- nrow(vc)
      s2 <- opt$value ^ 2 / (n - npar)
      std.errors <- sqrt(diag(vc) * s2)
      unc.alg <- "Standard errors computed from Hessian."
      
      # compute standard errors with backup algorithm
      if (all(is.na(std.errors))) {
        prior <- rbind(opt.df[1:5, -(1:2)], opt$par)
        v <- apply(prior, 2, var)
        std.errors <- sqrt(v * s2)
        message("Standard errors computed with backup algorithm.")
        unc.alg <-
          "Standard errors computed from variance of best optimization trials"
      }
      names(std.errors) <- c("mn", "mx", "sos", "rsp", "eos", "rau")
    }
    #    if (opt$convergence != 0) {
    #        opt$par[] <- NA
    #        xpred <- rep(NA, length(tout))
    #    } else {
    xpred <- .doubleLog(opt$par, tout)
    #    }
    xpred.out <- zoo(xpred, order.by = tout)
    if (plot) {
      llply(opt.l, function(opt) {
        xpred <- .doubleLog(opt$par, t)
        lines(t, xpred, col = "cyan")
      })
      lines(t, xpred, col = "blue", lwd = 2)
    }
    names(opt$par) <- c("mn", "mx", "sos", "rsp", "eos", "rau")
    fit.formula <- expression(mn + (mx - mn) * (1 / (1 + exp(-rsp *
                                                               (
                                                                 t - sos
                                                               ))) + 1 / (1 + exp(rau * (
                                                                 t - eos
                                                               ))) - 1))
    output <-
      list(predicted = xpred.out,
           params = opt$par,
           formula = fit.formula)
    if (hessian)
      output <- list(
        predicted = xpred.out,
        params = opt$par,
        formula = fit.formula,
        stdError = std.errors
      )
    
    return(output)
    ### The function returns a list with fitted values, parameters, fitting formula and standard errors if hessian is TRUE
  }

ex = function() {
  # Select one year of data.
  x <-
    as.vector(window(ndvi, start = c(1994, 1), end = c(1994, 12)))
  plot(x)
  
  # Fit double-logistic function to one year of data.
  fit <- FitDoubleLogBeck(x)
  lines(fit$predicted, col = "blue")
}
