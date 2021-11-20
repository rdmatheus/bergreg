#' Hypothesis Test for Dispersion Index in the BerG Regression (Under Development)
#'
#' \code{disp_test} performs a hypothesis test on the coefficients associated
#'     with the dispersion index in the BerG regression.
#'
#' @param object a 'bergreg' object.
#' @param cols an integer vector specifying which regressors must be tested.
#'    If \code{NULL} (default), a constant dispersion test is performed, that is,
#'    all coefficients, with the exception of the intercept, are tested to be
#'    equal to zero.
#' @param gamma0 vector of the null hypothesis; if \code{NULL} (default) a
#'    vector of zeros is passed to the constant dispersion test.
#' @param res_start initial guess for the maximization of the log-likelihood
#'  function restricted to the null hypothesis.
#' @param control optimization control parameters.
#' @param ... additional optimization parameters passed to the function \code{berg_control}.
#'
#' @return It returns a matrix with score, Wald, likelihood ratio, and gradient
#'  statistics with the respective p-values.
#' @export
#'
#' @examples
#' \dontrun{
#' # Loading the 'grazing' dataset
#' data(grazing)
#' head(grazing)
#'
#' # Fitting a double model
#' # birds ~ BerG(mu_i, phi_i)
#' # log(mu_i) = beta_0 + beta_1 when_i + beta_2 grazed_i
#' # log(mu_i) = gamma_0 + gamma_1 when_i + gamma_2 grazed_i
#'
#' fit <- bergreg(birds ~ when + grazed | when + grazed, grazing)
#'
#' # Test for constant dispersion
#' # H0: gamma_1 = gamma_2 = 0
#' disp_test(fit)
#' }
#'
disp_test <- function(object, cols = NULL, gamma0 = NULL,
                      res_start = NULL,
                      control = berg_control(...), ...){

  ## extract fitted information
  y <- if(is.null(object$y)){
    stats::model.response(stats::model.frame(object))
  } else{
    object$y
  }

  X <- stats::model.matrix(object$terms$mean, stats::model.frame(object))
  Z <- stats::model.matrix(object$terms$dispersion, stats::model.frame(object))
  p <- NCOL(X); k <- NCOL(Z);
  n <- length(y)

  link <- object$link$mean
  link.phi <- object$link$dispersion

  if(is.null(cols)){
    cols <- 2:k
  }

  q <- length(cols)

  if (is.null(gamma0)){
    gamma0 <- as.matrix(rep(0, q))
  }else{
    gamma0 <- as.matrix(gamma0)
  }

  # Partition in Z dispersion index model matrix
  Z1 <- as.matrix(Z[, -cols]); Z2 = as.matrix(Z[, cols])
  Z_aux <- cbind(Z1, Z2)

  # Unrestricted mle
  theta_hat <- mle_berg(y, X, Z_aux, link, link.phi, control)$est

  # Restricted mle
  eq_constraint <- function(theta){
    return(cbind(matrix(0, q, p + k - q), diag(rep(1, q)))%*%theta -
             gamma0)
  }

  eq_constraint_jac <- function(theta){
    return(cbind(matrix(0, q, p + k - q), diag(rep(1, q))))
  }

  if(!is.null(res_start)){
    control$start <- res_start
  }

  theta_tilde <- mle_berg(y, X, Z_aux, link, link.phi, control,
                          eq_constraint = eq_constraint,
                          eq_constraint_jac = eq_constraint_jac)$est


  #Ug2 <- U_berg(theta_tilde, y, X, Z_aux, link, link.phi)[(p + k - q + 1):(p + k)]
  #Kg2 <- solve(K_berg(theta_tilde, X, Z_aux, link, link.phi))[(p + k - q + 1):(p + k), (p + k - q + 1):(p + k)]

  # Tests statistics

  # Score
  S <- t(U_berg(theta_tilde, y, X, Z_aux, link, link.phi))%*%
    solve(K_berg(theta_tilde, X, Z_aux, link, link.phi))%*%
    U_berg(theta_tilde, y, X, Z_aux, link, link.phi)

  # Wald
  W <- t(theta_hat[(p + k - q + 1):(p + k)] - gamma0)%*%solve(
    solve(K_berg(theta_hat, X, Z_aux, link, link.phi))[(p + k - q + 1):(p + k), (p + k - q + 1):(p + k)]
  )%*%(theta_hat[(p + k - q + 1):(p + k)] - gamma0)

  # Likelihood ratio
  LR <- 2 * (ll_berg(theta_hat, y, X, Z_aux, link, link.phi) -
               ll_berg(theta_tilde, y, X, Z_aux, link, link.phi))

  # Gradient
  G <- t(U_berg(theta_tilde, y, X, Z_aux, link, link.phi))%*%
    (theta_hat - theta_tilde)

  #t(Ug2)%*%(theta_hat[(p + k - q + 1):(p + k)] - gamma0)

  out <- matrix(c(S, stats::pchisq(S, q, lower.tail = FALSE),
                  W, stats::pchisq(W, q, lower.tail = FALSE),
                  LR, stats::pchisq(LR, q, lower.tail = FALSE),
                  G, stats::pchisq(G, q, lower.tail = FALSE)), nrow = 2)
  rownames(out) <- c("value", "pvalue")
  colnames(out) <- c("Score","Wald","Lik. Ratio","Gradient")
  out
}
