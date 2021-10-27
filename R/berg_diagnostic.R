berg_envel <- function(object, type = c("quantile", "pearson", "response"),
                     nsim = 99, level = 0.95)
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
  type <- match.arg(type, c("quantile", "pearson", "response"))
  res <- stats::residuals(object, type = type)

  ## Envelope
  rs <- matrix(0, n, nsim)
  for(i in 1:nsim){

    # Simulated sample
    y.tilde <- rberg(n, mu, phi)

    # Estimatives
    #ctrl <- object$control
    ctrl <- berg_control()
    ctrl$hessian <- FALSE
    ctrl$start <- c(object$coefficients$mean,
                    object$coefficients$dispersion)

    fit <- suppressWarnings(mle_berg(y.tilde, X, Z,
                                   link = link,
                                   link.phi = link.phi,
                                   control = ctrl))

    ## Coefficients
    beta <- fit$est[1:p]
    gamma <- fit$est[1:k + p]

    ## Fitted values
    fit$fitted.values <- g(link)$inv(X%*%beta)
    fit$phi <- g(link.phi)$inv(Z%*%gamma)
    fit$residuals <- y.tilde - fit$fitted.values
    fit$y <- y.tilde

    # Empirical residuals
    rs[, i] <- sort(stats::residuals(fit, type = type))
  }

  alpha <- (1 - level)/2

  # alpha and 1 - alpha quantiles
  lower <- apply(rs, 1, stats::quantile, probs = alpha)
  upper <- apply(rs, 1, stats::quantile, probs = 1 - alpha)

  # Median
  md <- apply(rs, 1, stats::quantile, probs = 0.5)

  # Theoretical quantiles for the normal distribution
  qq <- stats::qqnorm(1:n, plot.it = FALSE)$x

  gp <- cbind(qq, sort(res), md, lower, upper)

  types <- c("quantile", "pearson", "response")
  Types <- c("Quantile residuals", "Pearson residuals", "Raw response residuals")
  Type <- Types[type == types]
  graphics::plot(gp[,1], gp[,2],
                 xlab = "Normal quantiles", ylab = Type, type = "n")
  graphics::polygon (c(gp[,1], sort(gp[,1], decreasing = T)),
                     c(gp[,4], sort(gp[,5], decreasing = T)),
                     col = "lightgray", border=NA)
  graphics::lines(gp[,1], gp[,3], lty = 2)
  graphics::points(gp[,1], gp[,2], pch = "+")
  graphics::box()
}
