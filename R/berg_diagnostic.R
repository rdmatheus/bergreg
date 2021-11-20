#' @name berg_envel
#'
#' @title  Envelope Plot for the BerG Regression Model
#'
#' @description   Provides the simulated envelope plot of the residuals (selectable
#' by \code{type}) resulting from the BerG regression fit.
#'
#' @param object a \code{'bergrm'} object.
#' @param x an \code{'envel_berg'} object.
#' @param type character; specifies which residual should be produced in the
#'     envelope plot. The available options are \code{"quantile"} (default),
#'     \code{"pearson"}, and \code{"response"} (raw residuals, y - mu).
#' @param nsim the number of replicates. The default is \code{nsim = 99}.
#' @param level level of the sample quantile of the residual used in the
#'     construction of confidence bands.
#' @param progressBar logical; if \code{TRUE}, a progress bar is displayed giving
#'     the progress of making the graph. It can slow down the function
#'     considerably in applications with a large sample size.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{berg_envel} returns an \code{"berg_envel"} object wichi consists of
#'     a list with the following components:
#' \describe{
#'   \item{residuals}{a list with the quantile, pearson, and raw residuals resulting from the
#'       fit of the BerG regression model.}
#'   \item{simulation}{a list whose components are matrices containing the
#'       ordered quantile, pearson, and raw residuals of the simulation for the
#'       plot envelope.}
#'
#'  The method \code{plot} makes the envelope plot.
#'  }
#'
NULL

#' @rdname berg_envel
#' @export
berg_envel <- function(object, nsim = 99, progressBar = TRUE, ...)
{
  ## Model specifications
  y <- if(is.null(object$y)) stats::model.response(stats::model.frame(object)) else object$y
  X <- if(is.null(object$x)) stats::model.matrix(object$terms$mean, stats::model.frame(object)) else object$y
  Z <- if(is.null(object$x)) stats::model.matrix(object$terms$dispersion, stats::model.frame(object)) else object$y
  p <- NCOL(X)
  k <- NCOL(Z)
  n <- length(y)

  link <- object$link$mean
  link.phi <- object$link$dispersion

  ## Fitted means and dispersions
  mu <- object$fitted.values
  phi <- object$phi

  ## Residuals
  rq <- stats::residuals(object, type = "quantile")
  rp <- stats::residuals(object, type = "pearson")
  rr <- stats::residuals(object, type = "response")

  ## Simulation
  rq_s <- matrix(0, n, nsim)
  rp_s <- matrix(0, n, nsim)
  rr_s <- matrix(0, n, nsim)

  pb <- utils::txtProgressBar(1, nsim, style = 3)
  for(i in 1:nsim){

    # Simulated sample
    y.tilde <- rberg(n, mu, phi)

    # Estimatives
    #ctrl <- object$control
    ctrl <- berg_control()
    ctrl$hessian <- FALSE
    ctrl$start <- c(object$coefficients$mean,
                    object$coefficients$dispersion)

    out <- suppressWarnings(mle_berg(y.tilde, X, Z,
                                     link = link,
                                     link.phi = link.phi,
                                     control = ctrl))

    ## Coefficients
    beta <- out$est[1:p]
    gamma <- out$est[1:k + p]

    ## Fitted values
    out$fitted.values <- g(link)$inv(X%*%beta)
    out$phi <- g(link.phi)$inv(Z%*%gamma)
    out$residuals <- y.tilde - out$fitted.values
    out$y <- y.tilde
    class(out) <- "bergreg"

    # Empirical residuals
    rq_s[, i] <- sort(stats::residuals(out, type = "quantile"))
    rp_s[, i] <- sort(stats::residuals(out, type = "pearson"))
    rr_s[, i] <- sort(stats::residuals(out, type = "response"))

    if(progressBar){ utils::setTxtProgressBar(pb, i)}
  }

  out <- list(residuals = list(quantile = rq,
                               pearson = rp,
                               response = rr),
              simulation = list(quantile = rq_s,
                                pearson = rp_s,
                                response = rr_s))

  class(out) <- "berg_envel"
  out

}

#' @rdname berg_envel
#' @export
print.berg_envel <- function(x, ...){
  cat("\nA 'berg_envel' object\n\n")
  utils::str(x)
}

#' @rdname berg_envel
#' @export
plot.berg_envel <- function(x, type = c("quantile", "pearson", "response"),
                            level = 0.95, ...){


  alpha <- (1 - level)/2
  type <- match.arg(type)

  if(type == "quantile"){
    res <- x$residuals$quantile
    n <- length(res)
    rs <- x$simulation$quantile

    # alpha and 1 - alpha quantiles
    lower <- apply(rs, 1, stats::quantile, probs = alpha)
    upper <- apply(rs, 1, stats::quantile, probs = 1 - alpha)

    # Median
    md <- apply(rs, 1, stats::quantile, probs = 0.5)

    # Theoretical quantiles for the normal distribution
    qq <- stats::qqnorm(1:n, plot.it = FALSE)$x

    gp <- cbind(qq, sort(res), md, lower, upper)
    graphics::plot(gp[,1], gp[,2],
                   xlab = "Normal quantiles", ylab = "Randomized quantile residuals", type = "n", ...)
    graphics::polygon (c(gp[,1], sort(gp[,1], decreasing = T)),
                       c(gp[,4], sort(gp[,5], decreasing = T)),
                       col = "lightgray", border=NA)
    graphics::lines(gp[,1], gp[,3], lty = 2)
    graphics::points(gp[,1], gp[,2], pch = "+")
    graphics::box()
  }

  if(type == "pearson"){

    res <- x$residuals$pearson
    n <- length(res)
    rs <- x$simulation$pearson

    # alpha and 1 - alpha quantiles
    lower <- apply(rs, 1, stats::quantile, probs = alpha)
    upper <- apply(rs, 1, stats::quantile, probs = 1 - alpha)

    # Median
    md <- apply(rs, 1, stats::quantile, probs = 0.5)

    # Theoretical quantiles for the normal distribution
    qq <- stats::qqnorm(1:n, plot.it = FALSE)$x

    gp <- cbind(qq, sort(res), md, lower, upper)
    graphics::plot(gp[,1], gp[,2],
                   xlab = "Normal quantiles", ylab = "Pearson residuals", type = "n", ...)
    graphics::polygon (c(gp[,1], sort(gp[,1], decreasing = T)),
                       c(gp[,4], sort(gp[,5], decreasing = T)),
                       col = "lightgray", border=NA)
    graphics::lines(gp[,1], gp[,3], lty = 2)
    graphics::points(gp[,1], gp[,2], pch = "+")
    graphics::box()
  }

  if(type == "response"){
    res <- x$residuals$response
    n <- length(res)
    rs <- x$simulation$response

    # alpha and 1 - alpha quantiles
    lower <- apply(rs, 1, stats::quantile, probs = alpha)
    upper <- apply(rs, 1, stats::quantile, probs = 1 - alpha)

    # Median
    md <- apply(rs, 1, stats::quantile, probs = 0.5)

    # Theoretical quantiles for the normal distribution
    qq <- stats::qqnorm(1:n, plot.it = FALSE)$x

    gp <- cbind(qq, sort(res), md, lower, upper)
    graphics::plot(gp[,1], gp[,2],
                   xlab = "Normal quantiles", ylab = "Raw residuals", type = "n", ...)
    graphics::polygon (c(gp[,1], sort(gp[,1], decreasing = T)),
                       c(gp[,4], sort(gp[,5], decreasing = T)),
                       col = "lightgray", border=NA)
    graphics::lines(gp[,1], gp[,3], lty = 2)
    graphics::points(gp[,1], gp[,2], pch = "+")
    graphics::box()
  }

  invisible(x)

}








