###################################################################
# Control optimization function in nloptr                         #
##################################################################
berg_control <- function(start = NULL,
                         start2 = NULL,
                         constant = 1e-7,
                         error = 1e-4,
                         algorithm = "NLOPT_LD_SLSQP",
                         method = "BFGS",...){

  rval <- list(start = start,
               start2 = start,
               constant = constant,
               error = error,
               algorithm = algorithm,
               method = method)

  rval <- c(rval, list(...))
  return(rval)
}

mle_berg <- function(y, X, Z, link = "log", link.phi = "log",
                     control = berg_control(...), eq_constraint = NULL,
                     eq_constraint_jac = NULL,
                     optimizer = c("nloptr", "optim"), ...){

  # Control list
  start     <- control$start
  constant  <- control$constant
  error     <- control$error
  algorithm <- control$algorithm
  method <- control$method

  # Link funtions

  # Dispersion link function
  g1 <- g(link)$fun
  g1.inv <- g(link)$inv
  g1. <- g(link)$deriv.

  # Dispersion link function
  g2 <- g(link.phi)$fun
  g2.inv <- g(link.phi)$inv
  g2. <- g(link.phi)$deriv.

  # Design matrices and necessary quantities
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  p <- NCOL(X)
  k <- NCOL(Z)
  n <- as.numeric(length(y))

  if(is.null(start)){
    bs <- solve(t(X)%*%X)%*%t(X)%*%g1(y + 0.1)
    mu. <- g1.inv(X%*%bs)
    sigma. <- stats::var(g1(y + 0.1))/(g1.(mu.)^2)
    phi. <- sigma. / mu. + abs(mu. - 1)
    gs <- solve(t(Z)%*%Z)%*%t(Z)%*%g2(phi.)
    start <- c(bs, gs)
  }

  # Log-likelihood
  ll <- function(par) -ll_berg(par, y, X, Z, link, link.phi)

  # Score function
  U <- function(par) -U_berg(par, y, X, Z, link, link.phi)

  # Constraints and its Jacobian
  hj <- function(par){


    # Relations
    mu  <- g1.inv(X%*%par[1:p])
    phi <- g2.inv(Z%*%par[(p + 1):(p + k)])

    # Function h (h (theta) < 0 if theta is admissible)
    hj <- rbind(- phi + mu - 1 + constant, - phi - mu + 1 + constant)

    return(hj)

  }

  Jh <- function(par){

    # Relations
    mu  <- g1.inv(X%*%par[1:p])
    phi <- g2.inv(Z%*%par[(p + 1):(p + k)])

    # Diagonal matrix
    D1 <- diag(as.numeric(1 / g1.(mu)))
    D2 <- diag(as.numeric(1 / g2.(phi)))

    # Jacobian matrix
    Jh <- matrix(NA, 2 * n, p + k)
    Jh[1:n, 1:p] <- D1%*%X
    Jh[1:n, (p + 1):(p + k)] <- - D2%*%Z
    Jh[(n + 1):(2 * n), 1:p] <- - D1%*%X
    Jh[(n + 1):(2 * n),(p + 1):(p + k)] <- - D2%*%Z


    return(Jh)
  }

  optimizer <- match.arg(optimizer, c("nloptr", "optim"))

  if (optimizer == "nloptr"){

    if (!is.null(eq_constraint)){
      eq_constraint <- match.fun(eq_constraint)
      eq_constraint_jac <- match.fun(eq_constraint_jac)
    }

    est <- nloptr::nloptr(x0 = start,
                          eval_f = ll,
                          eval_grad_f = U,
                          eval_g_ineq = hj,
                          eval_jac_g_ineq = Jh,
                          eval_g_eq = eq_constraint,
                          eval_jac_g_eq = eq_constraint_jac,
                          opts = list("algorithm" = algorithm,
                                      "xtol_rel"= error))$solution

    logLik <- -ll(est)

    mu.h  <- g1.inv(X%*%est[1:p])
    phi.h <- g2.inv(Z%*%est[(p + 1):(p + k)])

    feasible <- all(phi.h > abs(mu.h - 1))
  }

  if (optimizer == "optim"){

    control$start <- control$start2 <- control$constant <-
      control$error <- control$optimizer <- control$algorithm <-
      control$method <- NULL

    est <- suppressWarnings(stats::optim(par = start,
                                         fn = ll,
                                         gr = U,
                                         method = method,
                                         control = control)$par)
    logLik <- -ll(est)

    mu.h  <- g1.inv(X%*%est[1:p])
    phi.h <- g2.inv(Z%*%est[(p + 1):(p + k)])

    feasible <- all(phi.h > abs(mu.h - 1))
  }

  return(list(est = est, logLik = logLik, feasible = feasible))
}
