#' @name bergreg-methods
#' @title Methods for 'bergreg' objects
#' @param x,object an object of class \code{bergreg}.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
NULL

# Print
#' @rdname bergreg-methods
#' @export
print.bergreg <- function(x, ...)
{

  #X <- stats::model.matrix(x$terms$mean, stats::model.frame(x))
  #Z <- stats::model.matrix(x$terms$dispersion, stats::model.frame(x))
  p <- length(x$coefficients$mean)
  k <- length(x$coefficients$dispersion)
  n <- x$nobs

  cat("\nCall:\n")
  print(x$call)
  #if(!x$converged) {
  #  cat("\nmodel did not converge\n")
  #}else{
    cat("\nmean Coefficients:\n")
    print((x$coefficients)$mean)
    cat("\ndispersion Coefficients:\n")
    print((x$coefficients)$dispersion)
    cat("\nIn addition, Log-lik value:",x$logLik,
        "\nAIC:", 2 * (p + k - x$logLik),
        "and BIC:", log(n) * (p + k) - 2 * x$logLik)
  #}

  invisible(x)
}

# Summary
#' @rdname bergreg-methods
#' @export
summary.bergreg <- function(object, ...)
{
  X <- stats::model.matrix(object$terms$mean, stats::model.frame(object))
  Z <- stats::model.matrix(object$terms$dispersion, stats::model.frame(object))
  p <- ncol(X)
  k <- ncol(Z)
  n <- object$nobs

  ## Links
  link <- object$link$mean
  link.phi <- object$link$dispersion

  ## Names
  if (p > 1) {
    mean_names <- colnames(X)
  }else{
    if (link == "identity"){
      mean_names <- "mu"
    }else{
      mean_names <- "g1(mu)"
    }
  }

  if (k > 1) {
    disp_names <- colnames(Z)
  }else{
    if (link.phi == "identity"){
      disp_names <- "phi"
    }else{
      disp_names <- "g2(phi)"
    }
  }

  ## Summary for quantile residuals
  res <- stats::residuals(object, type = "quantile")
  skewness <- mean((res - mean(res))^3) / (stats::sd(res)^3)
  kurtosis <- mean((res - mean(res))^4) / (stats::sd(res)^4)
  TAB.residuals <- round(cbind(mean(res), stats::sd(res),
                               skewness, kurtosis), 6)
  colnames(TAB.residuals) <- c("Mean", "Std. dev.", "Skewness", "Kurtosis")
  rownames(TAB.residuals) <- " "

  # Summary for mu
  est.beta <- object$coef$mean
  se.beta <- sqrt(diag(object$vcov))[1:p]
  zval.beta <- est.beta/se.beta
  pval.beta <- 2 * stats::pnorm(abs(zval.beta), lower.tail = FALSE)

  TAB.mu <- cbind(Estimate = est.beta,
                  `Std. Error` = se.beta,
                  `z value` = zval.beta,
                  `Pr(>|z|)` = pval.beta)
  rownames(TAB.mu) <- mean_names

  # Summary for sigma
  est.gamma <- object$coef$disp
  se.gamma <- sqrt(diag(object$vcov)[1:k + p])
  zval.gamma <- est.gamma/se.gamma
  pval.gamma <- 2*stats::pnorm(abs(zval.gamma), lower.tail = FALSE)

  TAB.phi <- cbind(Estimate = est.gamma,
                   `Std. Error` = se.gamma,
                   `z value` = zval.gamma,
                   `Pr(>|z|)` = pval.gamma)
  rownames(TAB.phi) <- disp_names

  # Links
  link <- object$link$mean
  link.phi <- object$link.sigma

  # logLik, AIC, and BIC
  logLik <- object$logLik
  AIC <- 2 * (p + k - logLik)
  BIC <- log(n) * (p + k) - 2 * logLik

  out <- list(call = object$call,
              residuals = TAB.residuals,
              link = object$link$mean,
              link.phi = object$link$dispersion,
              mean = TAB.mu,
              dispersion = TAB.phi,
              logLik = logLik,
              AIC = AIC,
              BIC = BIC)

  class(out) <- "summary.bergreg"
  out
}

# Print summary
#' @rdname bergreg-methods
#' @export
print.summary.bergreg <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nSummary for quantile residuals:\n")
  print(x$residuals)
  cat("\n---------------------------------------------------------------\n")
  cat("Mean model with", x$link, "link:")
  cat("\n---------------------------------------------------------------\n")
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$mean, signif.legend = FALSE)
  cat("\n---------------------------------------------------------------\n")
  cat("Dispersion model with", x$link.phi, "link:")
  cat("\n---------------------------------------------------------------\n")
  #cat("\nLink function:",x$link.sigma,"\n")
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$dispersion)
  cat("\n---------------------------------------------------------------")
  cat("\nIn addition, Log-lik value:",x$logLik,
      "\nAIC:", x$AIC,
      "and BIC:", x$BIC)

  invisible(x)
}

# Parameter estimates
#' @rdname bergreg-methods
#' @export
#' @param model a character indicating which parameter coefficients are
#'   required, parameters for the \code{"mean"} or for the
#'   \code{"dispersion"} model. If \code{"full"} (default), a list with
#'   coefficients for the \code{mean} and for the \code{dispersion}
#'   model is returned.
coef.bergreg <- function(object,
                       model = c("full", "mean", "dispersion"), ...) {
  model <- match.arg(model)
  switch(model,
         "full"        = list(
           mean       = (object$coef)$mean,
           dispersion = (object$coef)$disp),
         "mean"       = (object$coef)$mean,
         "dispersion" = (object$coef)$disp)
}

#  Variance-covariance matrix
#' @rdname bergreg-methods
#' @export
vcov.bergreg <- function(object,
                       model = c("full", "mean", "dispersion"), ...) {

  model <- match.arg(model)
  covm <- object$vcov

  X <- stats::model.matrix(object$terms$mean, stats::model.frame(object))
  Z <- stats::model.matrix(object$terms$dispersion, stats::model.frame(object))
  p <- ncol(X)
  k <- ncol(Z)

  ## Names
  if (p > 1) {
    mean_names <- colnames(X)
  }else{
    if (object$link$mean == "identity"){
      mean_names <- "mu"
    }else{
      mean_names <- "g1(mu)"
    }
  }

  if (k > 1) {
    disp_names <- colnames(Z)
  }else{
    if (object$link$dispersion == "identity"){
      disp_names <- "phi"
    }else{
      disp_names <- "g2(phi)"
    }
  }

  switch(model,
         "mean" = {
           covm <- covm[seq.int(length.out = p),
                        seq.int(length.out = p), drop = FALSE]
           colnames(covm) <- rownames(covm) <- mean_names
           covm
         },
         "dispersion" = {
           covm <- covm[seq.int(length.out = k) + p,
                        seq.int(length.out = k) + p, drop = FALSE]
           colnames(covm) <- rownames(covm) <- disp_names
           covm
         },
         "full" = {
           colnames(covm) <- rownames(covm) <- c(mean_names,
                                                 disp_names)
           covm
         })

}


#' Predict Method for BerG Fits
#'
#' Obtains predictions from a fitted BerG regression object.
#'
#' @param object an \code{'bergreg'} object.
#' @param newdata optionally, a data frame in which to look for variables
#'     with which to predict. If omitted, the fitted linear predictors are
#'     used.
#' @param type the type of prediction required. The default is on the scale of
#'     the response variable \code{("response")}, that is, the fitted values
#'     (fitted means). The alternative \code{"link"} is on the scale of
#'     the linear predictors, that is, predictions are of log-means. The
#'     specification \code{"dispersion"} provides the fitted dispersion
#'     indixes, while \code{"variance"} provides the fitted variances. Finally,
#'     the option \code{"quantile"} gives the fitted quantiles in the order
#'     specified via \code{'at'}.
#' @param at the order of the quantile to be predicted if
#'     \code{type = "quantile"}. The default is to predict the median,
#'     that is, \code{at = 0.5}.
#' @param na.action function determining what should be done with missing
#'     values in newdata. The default is to predict \code{NA}.
#' @param ...  arguments passed to or from other methods.
#'
#' @return A vector of predictions.
#' @export
#'
predict.bergreg <- function(object, newdata = NULL,
                           type = c("response", "link", "dispersion",
                                    "variance", "quantile"),
                           at = 0.5,
                           na.action = stats::na.pass, ...)
{
  type <- match.arg(type)

  qfun <- function(at, mu, phi) {
    rval <- sapply(at, function(p) qberg(rep(p, length(mu)), mu, phi))
    if(length(at) > 1L) {
      if(NCOL(rval) == 1L)
        rval <- matrix(rval, ncol = length(at),
                       dimnames = list(unique(names(rval)), NULL))

      colnames(rval) <- paste("q_", at, sep = "")
    } else {
      rval <- drop(rval)
    }
    rval
  }


  if(missing(newdata)) {

    rval <- switch(type,
                   "response" = {
                     as.numeric(object$fitted.values)
                   },
                   "link" = {
                     as.numeric(g(object$link$mean)$fun(object$fitted.values))
                   },
                   "dispersion" = {
                     object$phi
                   },

                   "variance" = {
                     mu <- object$fitted.values
                     phi <- object$phi
                     as.numeric(phi * mu)
                   },
                   "quantile" = {
                     mu <- as.numeric(object$fitted.values)
                     phi <- as.numeric(object$phi)

                     qfun(at, mu, phi)
                   })

    return(rval)

  } else {

    tnam <- switch(type,
                   "response" = "mean",
                   "link" = "mean",
                   "dispersion" = "dispersion",
                   "variance" = "full",
                   "quantile" = "full")

    mf <- stats::model.frame(stats::delete.response(object$terms[[tnam]]),
                             newdata, na.action = na.action)
    newdata <- newdata[rownames(mf), , drop = FALSE]

    if(type %in% c("response", "link", "variance", "quantile"))
      X <- stats::model.matrix(stats::delete.response(object$terms$mean), mf)

    if(type %in% c("dispersion", "variance", "quantile"))
      Z <- stats::model.matrix(object$terms$dispersion, mf)

    rval <- switch(type,
                   "response" = {
                     as.numeric(g(object$link$mean)$inv(
                       drop(X %*% object$coefficients$mean)))
                   },
                   "link" = {
                     mu <- as.numeric(g(object$link$mean)$inv(
                       drop(X %*% object$coefficients$mean)))

                     g(object$link$mean)$fun(mu)
                   },
                   "dispersion" = {
                     as.numeric(g(object$link$dispersion)$inv(
                       drop(Z %*% object$coefficients$dispersion)))
                   },
                   "variance" = {
                     mu <- g(object$link$mean)$inv(
                       drop(X %*% object$coefficients$mean))

                     phi <- g(object$link$dispersion)$inv(
                       drop(Z %*% object$coefficients$dispersion))
                     as.numeric(mu * phi)

                   },
                   "quantile" = {
                     mu <- g(object$link$mean)$inv(
                       drop(X %*% object$coefficients$mean))

                     phi <- g(object$link$dispersion)$inv(
                       drop(Z %*% object$coefficients$dispersion))

                     qfun(at, mu, phi)
                   }
    )
    return(rval)

  }
}


# Log-likelihood
#' @rdname bergreg-methods
#' @export
logLik.bergreg <- function(object, ...) {
  structure(object$logLik,
            df = object$nobs - object$df.residual,
            class = "logLik")
}


# AIC
#' @export
#' @rdname bergreg-methods
AIC.bergreg <- function(object, ..., k = 2) {

  AIC <- - 2 * object$logLik + k *
    (length(object$coefficients$mean) + length(object$coefficients$dispersion))

  class(AIC) <- "AIC"
  return(AIC)
}

# Residuals
#' @name residuals.bergreg
#' @title Extract Model Residuals for a BerG Regression
#'
#' @param object an \code{'bergreg'} object.
#' @param type character; specifies which residual should be extracted.
#'     The available arguments are "quantile" (default), "pearson",
#'     and "response" (raw residuals, y - mu).
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
residuals.bergreg <- function(object,
                            type = c("quantile", "pearson", "response"), ...)
{
  ## raw response residuals and desired type
  res <- object$residuals
  type <- match.arg(type, c("quantile", "pearson", "response"))
  if(type == "response") return(res)

  ## extract fitted information
  y <- if(is.null(object$y)){
         stats::model.response(stats::model.frame(object))
       } else{
         object$y
       }

  mu <- stats::fitted(object)
  phi <- object$phi

  rqr_berg <- function(y){
    n <- length(y)
    a <- vector()
    b <- vector()
    u <- vector()
    for(i in 1:n){
      a[i] <- pberg(y[i] - 1, mu[i], phi[i])
      b[i] <- pberg(y[i], mu[i], phi[i])
      u[i] <- stats::runif(1, a[i], b[i])
    }
    return(stats::qnorm(u))
  }

  res <- switch(type,
                "pearson" = {
                  as.numeric(res / sqrt(phi * mu))
                },

                #"deviance" = {
                #  ll <- function(mu, sigma) dugamma(y, mu, sigma, log.p = TRUE)
                #  sign(res) * sqrt(2 * abs(ll(y, sigma) - ll(mu, sigma)))
                #},

                "quantile" = {
                  rqr_berg(y)
                })

  res
}

## Model frame
model.frame.bergreg <- function(formula, ...) {
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

## Model matrix
#' @export
#' @rdname bergreg-methods
model.matrix.bergreg <- function(object,
                               model = c("mean", "dispersion"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
  else stats::model.matrix(object$terms[[model]], stats::model.frame(object))
  return(rval)
}

# Plot
#' Diagnostic plots for the BerG Regression
#'
#' Five plots (selectable by \code{which}) are currently available:
#' a plot of residuals against fitted values, a plot of residuals against
#' the indexes, a Normal Q-Q plot, a barplot with comparisons of the
#' observed and fitted frequencies, and a plot of the sample autocorelations
#' of the residuals.
#'
#' @param x an object of class \code{bergreg}.
#' @param type character; specifies which residual should be produced
#'     in the summary plot. The available arguments are "quantile", "pearson",
#'     and "response" (raw residuals, y - mu).
#' @param which numeric; if a subset of the plots is required, specify a subset
#'     of the numbers 1:5.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
plot.bergreg <- function(x, type = c("quantile", "pearson", "response"),
                         which = 1:2,
                         ask = prod(graphics::par("mfcol")) < length(which) &&
                         grDevices::dev.interactive(),
                         ...)
{
  if(!is.numeric(which) || any(which < 1) || any(which > 5))
    stop("`which' must be in 1:5")

  ## Reading
  type <- match.arg(type)
  res <- stats::residuals(x, type = type)

  ## Legends
  types <- c("quantile", "pearson", "response")
  Types <- c("Quantile residuals", "Pearson residuals", "Raw response residuals")
  Type <- Types[type == types]

  ## Graphical parameters setting
  if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op))
  }

  ## Plots to shown
  show <- rep(FALSE, 5)
  show[which] <- TRUE

  ## Residuals versus Fitted values
  if (show[1]){
    #y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
    graphics::plot(stats::fitted(x), res,
                   xlab = "Fitted values", ylab = Type, pch = "+")
    graphics::abline(h = 0, col = "gray50", lty = 3)
  }

  ## Residuals versus index observation
  if (show[2]){
    n <- x$nobs
    graphics::plot(1:n, res, xlab = "Index", ylab = Type, pch = "+")
    graphics::abline(h = 0, col= "gray50", lty = 3)
  }

  ## Normal probability plot
  if(show[3]) {
    stats::qqnorm(res, main = "Normal Q-Q Plot",
           xlab = "Theoretical Quantiles",
           ylab = "Sample Quantiles",
           plot.it = TRUE,
           frame.plot = TRUE, pch =  "+")
    graphics::abline(0, 1, col = "gray50", lty = 2)
  }

  ## Expected frequencies
  if(show[4]) {
    y <- if(is.null(x$y)) stats::model.response(stats::model.frame(x)) else x$y
    obs <- as.numeric(table(y))
    xcoord <- graphics::barplot(table(y), xlab = "Count", ylab = "Frequency",
                   ylim = c(0, max(max(obs), max(x$freq)) + 0.5))
    graphics::points(xcoord, x$freq, col = 2, type = "b", pch = 16)
    graphics::legend("topright", "Fitted frequencies", col = 2, lty = 1, pch = 16,
           bty = "n")
  }

  ## ACF of residuals
  if(show[5]) {
    stats::acf(res, main = " ", xlab = "Lags",
               ylab = paste("Sample ACF of", type, "residuals"))
  }

  invisible(x)

}
