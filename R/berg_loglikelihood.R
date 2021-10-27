###################################################################
# Link function                                                   #
##################################################################
g <- function(link){

  switch(link,

         identity = {
           fun <- function(theta) theta
           inv <- function(eta) eta
           deriv. <- function(theta) rep.int(1, length(theta))
           deriv.. <- function(theta) rep.int(0, length(theta))
           valideta <- function(eta) TRUE
         },

         log = {
           fun <- function(theta) log(theta)
           inv <- function(eta) pmax(exp(eta), .Machine$double.eps)
           deriv. <- function(theta) 1 / theta
           deriv.. <- function(theta) -1 / (theta^2)
           valideta <- function(eta) TRUE
         },

         sqrt = {
           fun <- function(theta) sqrt(theta)
           inv <- function(eta) eta^2
           deriv. <- function(theta) 1 / (2 * sqrt(theta))
           deriv.. <- function(theta) -1 / (4 * (theta ^ (3 / 2)))
           valideta <- function(eta) all(is.finite(eta)) && all(eta > 0)
         },

         stop(gettextf("link %s not available", sQuote(link)), domain = NA))

  environment(fun) <- environment(inv) <- environment(deriv.) <-
  environment(deriv..) <- environment(valideta) <- asNamespace("stats")

  structure(list(fun = fun, inv = inv, deriv. = deriv.,
                 deriv.. = deriv.., valideta = valideta,
                 name = link))
}


###################################################################
# Log-lokelihood function                                         #
##################################################################
ll_berg <- function(par, y, X, Z, link = "log", link.phi = "log"){

  # Link functions
  g1.inv <- g(link)$inv
  g2.inv <- g(link.phi)$inv

  # Design matrices and necessary quantities
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- NCOL(X); k <- NCOL(Z)

  # Relations
  mu <- g1.inv(X%*%par[1:p])
  phi <- g2.inv(Z%*%par[1:k + p])

  if(any(!is.finite(mu)) | any(!is.finite(phi)) | any(phi <= abs(mu - 1)) | any(mu < 0)){
    NaN
  }else {
    sum(dberg(y, mu, phi, log.p = TRUE))
  }
}


###################################################################
# Score function                                                  #
##################################################################
U_berg <- function(par, y, X, Z, link = "log", link.phi = "log"){

  # Link functions
  g1.inv <- g(link)$inv
  g1. <- g(link)$deriv.
  g1.. <- g(link)$deriv..

  g2.inv <- g(link.phi)$inv
  g2. <- g(link.phi)$deriv.
  g2.. <- g(link.phi)$deriv..

  # Design matrices and necessary quantities
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- NCOL(X); k <- NCOL(Z)

  # Relations
  mu <- g1.inv(X%*%par[1 : p])
  phi <- g2.inv(Z%*%par[(p + 1):(p + k)])

  id1 <- which(y == 0 & mu > 0 & phi > abs(mu - 1))
  id2 <- which(y > 0 & mu > 0 & phi > abs(mu - 1))

  u1 <- u2 <- vector("numeric", length = length(y))

  u1[id1] <- -2 * exp(log(1 + phi[id1]) - log(1 - mu[id1] + phi[id1]) - log(1 + mu[id1] + phi[id1]))
  u1[id2] <-  1 / mu[id2] + 2 * (y[id2] - mu[id2] - phi[id2])/
                                exp(log(mu[id2] + phi[id2] - 1) +
                                    log(mu[id2] + phi[id2] + 1))

  u2[id1] <- 2 * exp(log(mu[id1]) - log(1 - mu[id1] + phi[id1]) - log(1 + mu[id1] + phi[id1]))
  u2[id2] <- 2 * (y[id2] - mu[id2] - phi[id2]) /
                  exp(log(mu[id2] + phi[id2] - 1) +
                      log(mu[id2] + phi[id2] + 1))

  D1 <- diag(as.numeric(1/g1.(mu)))
  D2 <- diag(as.numeric(1/g2.(phi)))

  Ub <- t(X)%*%D1%*%u1; Ug = t(Z)%*%D2%*%u2

  c(Ub, Ug)
}

###################################################################
# Fisher information matrix                                       #
##################################################################
K_berg <- function(par, X, Z, link = "log", link.phi = "log", inverse = FALSE){

  # Link functions
  g1.inv <- g(link)$inv
  g1. <- g(link)$deriv.
  g1.. <- g(link)$deriv..

  g2.inv <- g(link.phi)$inv
  g2. <- g(link.phi)$deriv.
  g2.. <- g(link.phi)$deriv..

  # Design matrices and necessary quantities
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- NCOL(X); k <- NCOL(Z)

  # Relations
  mu <-  g1.inv(X%*%par[1 : p])
  phi <- g2.inv(Z%*%par[1:k + p])

  D1 <- diag(as.numeric(1 / g1.(mu)))
  D2 <- diag(as.numeric(1 / g2.(phi)))

  K <- matrix(rep(NA, (p + k) * (p + k)), p + k, p + k)

  W1 <- diag(as.numeric(

    2 * (2 * exp(log(mu) + log(1 + phi) - log(1 - mu + phi) - 3 * log(1 + mu + phi)) -
         2 * exp(log(mu) + log(1 + (mu + phi)^2) - 2 * log(mu + phi - 1) - 3 * log(mu + phi + 1)) +
         exp(log(mu) + log(mu + phi) - 2 * log(mu + phi - 1) - 2 * log(mu + phi + 1)) +
         exp(-log(mu) - log(mu + phi + 1))) / g1.(mu) +

    2 * (-exp(log(1 + phi) - 2 * log(1 + mu + phi)) + exp(-log(1 + mu + phi)) -
          2 * exp(log(mu) + log(mu + phi) - log(mu + phi - 1) - 2 * log(mu + phi + 1)) +
          exp(log(mu) - log(mu + phi - 1) - log(mu + phi + 1))) * (g1..(mu) / (g1.(mu)^2))

    ))

  W2 <- diag(as.numeric(

    -2 * exp(log(mu^2 + (1 + phi)^2) - log(1 - mu + phi) - 3 * log(1 + mu + phi)) +
     4 * exp(log(mu) + log(mu + phi) - 2 * log(mu + phi - 1) - 2 * log(mu + phi + 1)) -
     4 * exp(log(mu) + log(1 + (mu + phi)^2) - 2 * log(mu + phi - 1) - 3 * log(mu + phi + 1))

  ))

  W3 <- diag(as.numeric(

    4 * (exp(log(mu) + log(1 + phi) - log(1 - mu + phi) - 3 * log(1 + mu + phi)) +
         exp(log(mu) + log(mu + phi) - 2 * log(mu + phi - 1) - 2 * log(mu + phi + 1)) -
         exp(log(mu) + log(1 + (mu + phi)^2) - 2 * log(mu + phi - 1) - 3 * log(mu + phi + 1))) /
      g2.(phi) +

      2 * (exp(log(mu) - 2 * log(mu + phi + 1)) +
           exp(log(mu) - log(mu + phi - 1) - log(mu + phi + 1)) -
           2 * exp(log(mu) + log(mu + phi) - log(mu + phi - 1) - 2 * log(mu + phi + 1))) *
      (g2..(phi) / (g2.(phi)^2))

  ))

  Kbb <- t(X)%*%W1%*%D1%*%X
  Kbg <- t(X)%*%D1%*%W2%*%D2%*%Z
  Kgb <- t(Kbg)
  Kgg <- t(Z)%*%W3%*%D2%*%Z

  K[1:p,1:p] <- Kbb; K[1:p,(p + 1):(p + k)] <- Kbg; K[(p + 1):(p + k),(1:p)] <- Kgb
  K[(p + 1):(p + k),(p + 1):(p + k)] <- Kgg

  if(inverse) chol2inv(chol(K)) else K

}
