#' @name berg
#'
#' @title The BerG Distribution
#'
#' @description Probability mass function, distribution function,
#' quantile function, and random generation for the BerG distribution
#' with parameters \code{mu} and \code{phi}.
#'
#' @param x vector of non-negative integer quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu numeric; non-negative mean.
#' @param phi numeric; the dispersion index (greather than \code{abs(mu - 1)}),
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P(X <= x)}, otherwise, \code{P(X > x)}.
#' @param log.p 	logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#'
#' @return \code{dberg} returns the probability function, \code{pberg}
#' gives the distribution function, \code{qberg} gives the quantile function,
#' and \code{rberg} generates random observations.
#'
#' @details This set of functions represents the probability function, the cumulative distribution
#' function, quantile function, and a random number generator for the BerG distribution
#' parameterized in terms of the mean and the dispersion index. This new
#' parameterization was proposed by Bourguignon, M. & Medeiros, R. (2021).
#'
#' @references Bourguignon, M. & Weiss, C. (2017). An INAR(1) process for modeling count
#' time series with equidispersion, underdispersion and overdispersion. Test, 26, 847--868.
#'
#' @references Bourguignon, M. & Medeiros, R. (2021). A simple and useful regression model for fitting count data
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' # BerG observations as categorical data:
#' Ni <- rberg(50, mu = 4, phi = 3.2); table(factor(Ni, 0:max(Ni)))
#' plot(prop.table(table(Ni)))
#' points(sort(unique(Ni)), dberg(sort(unique(Ni)), mu = 4, phi = 3.2), pch=16, col="red")
#'
#' # Probability function:
#'
#' # Parameters
#' mu = c(2.3, 3, 3.9); phi = 3
#'
#' plot(0:17-0.25,dberg(0:17,phi,mu[1]),type="h",xlab=" ",ylab="Pmf",
#' col="gray70",lwd=5);mtext("y",side=1,line=2.5)
#' segments(0:17,rep(0,length(0:17)),0:17,dberg(0:17,phi,mu[2]),col="gray40",lwd=5)
#' segments(0:17+0.25,rep(0,length(0:17)),0:17+0.25,dberg(0:17,phi,mu[3]),col="gray20",lwd=5)
#'
#' legend(13,0.32,legend =c(expression(mu==2.3),
#'           expression(mu==3), expression(mu==3.9)),
#'           lty=1,lwd=5,col=c("gray70","gray40","gray20"),bty="n")
#'
NULL

#' @rdname berg
#' @export
dberg <- function(x, mu, phi, log.p = FALSE){

  if (is.vector(x))
    x <- matrix(x, ncol = length(x))

  n <- dim(x)[2]
  mu <- matrix(mu, ncol = n)
  phi <- matrix(phi, ncol = n)

  pmf <- matrix(-Inf, ncol = n)

  # NaN index
  NaNid <- which(mu <= 0 | phi <= abs(mu - 1), arr.ind = TRUE)
  pmf[NaNid] <- NaN

  # Positive pmf index
  id1 <- which(x == 0 & mu > 0 & phi > abs(mu - 1) & !is.nan(pmf), arr.ind = TRUE)
  id2 <- which(x > 0 & mu > 0 & phi > abs(mu - 1) & !is.nan(pmf), arr.ind = TRUE)

  pmf[id1] <- log(1 - mu[id1] + phi[id1]) - log(1 + mu[id1] + phi[id1])
  pmf[id2]  <- log(4) + log(mu[id2]) + (x[id2] - 1) * log(mu[id2] + phi[id2] - 1) -
                (x[id2] + 1) * log(mu[id2] + phi[id2] + 1)

  if(nrow(pmf) == 1L)
    pmf <- as.vector(pmf)

  if(log.p) pmf else exp(pmf)
}

#' @rdname berg
#' @export
pberg <- function(q, mu, phi, lower.tail = TRUE){

  q <- floor(q)

  if (is.vector(q))
    q <- matrix(q, ncol = length(q))

  n <- dim(q)[2]
  mu <- matrix(mu, ncol = n)
  phi <- matrix(phi, ncol = n)

  cdf <- matrix(0, ncol = n)

  # NaN index
  NaNid <- which(mu <= 0 | phi <= abs(mu - 1), arr.ind = TRUE)
  cdf[NaNid] <- NaN

  # Positive cdf index
  id <- which(q >= 0 & mu > 0 & phi > abs(mu - 1) & !is.nan(cdf), arr.ind = TRUE)
  cdf[id]  <- exp(log(1 - mu[id] + phi[id]) - log(1 + mu[id] + phi[id])) +
      exp(log(2 * mu[id]) - log(1 + mu[id] + phi[id]) +
    log(1 - ((mu[id] + phi[id] - 1) / (mu[id] + phi[id] + 1))^q[id]))

  if(nrow(cdf) == 1L)
    cdf <- as.vector(cdf)

  if(lower.tail) cdf else 1 - cdf
}

#' @rdname berg
#' @export
qberg <- function(p, mu, phi, lower.tail = TRUE){

  if(!lower.tail)
    p <- 1 - p

  if (is.vector(p))
    p <- matrix(p, ncol = length(p))

  n <- dim(p)[2]
  mu <- matrix(mu, ncol = n)
  phi <- matrix(phi, ncol = n)

  q <- matrix(0, ncol = n)

  # NaN index
  NaNid <- which(mu <= 0 | phi <= abs(mu - 1), arr.ind = TRUE)
  q[NaNid] <- NaN

  p0 <- exp(log(1 - mu + phi) - log(1 + mu + phi))

  id <- which(p >= p0 & mu > 0 & phi > abs(mu - 1) & !is.nan(q), arr.ind = TRUE)
  q[id] <- ceiling(
    (log(1 - p[id]) + log(1 + mu[id] + phi[id]) - log(2 * mu[id])) /
    (log(mu[id] + phi[id] - 1) - log(1 + mu[id] + phi[id]))
  )

  if(nrow(q) == 1L) as.vector(q) else q

}

#' @rdname berg
#' @export
rberg <- function(n, mu, phi){

  if ((any(mu<=0)) || (any(phi<=0)))
    stop("The parameters must be positives")
  if (any(phi < abs(mu - 1)))
    stop("Constraints are not satisfied")

  u <- stats::runif(n)
  return(qberg(u, mu, phi))
}

