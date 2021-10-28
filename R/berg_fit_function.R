#' @name bergreg
#'
#' @title BerG Regression for Count Data
#'
#' @description Fit of the BerG regression model via maximum
#' likelihood for a new parameterization of this distribution
#' that is indexed by the mean and the dispersion index.
#'
#' @param formula a simbolic description of the model, of type
#'  \code{y ~ x} for covariates in the mean submodel only or \code{y ~ x | z}
#'      to enter covariates in the dispersion index submodel. See details below.
#' @param data an optional data frame containing the variables in the formula.
#'     By default the variables are taken from environment(formula).
#' @param link,link.phi character; specification of the link function in the
#'     mean and dispersion index submodels. The links \code{"log"} (default),
#'     \code{"sqrt"}, and \code{"identity"} can be used.
#' @param disp.test logical; if \code{TRUE}, the function \code{bergreg} returns the
#'     test for constant dispersion.
#' @param optimizer character; specification of the optimization algorithm.
#'     By default, estimates are obtained using sequential quadratic programming
#'     method via nloptr package (\code{"nloptr"}). However, it is also possible to obtain the
#'     estimates via the \code{optim} function (\code{"optim"}).
#' @param y,x logicals; if \code{TRUE} the corresponding components of the fit,
#'      response and model matrices, are returned.
#' @param control a list of control arguments specified via \code{berg_control}
#'     (under development).
#' @param ... arguments passed to \code{berg_control} (under development).
#'
#' @return  The \code{bergreg} function returns an object of class "bergreg",
#'  which consists of a list with the following components:
#' \describe{
#'   \item{coefficients}{a list containing the elements "mean" and "dispersion,"
#'       which consists of the estimates of the coefficients associated with the
#'       mean and the dispersion index, respectively.}
#'   \item{fitted.values}{a vector with the fitted means.}
#'   \item{phi}{a vector with the fitted dispersion indexes.}
#'   \item{residuals}{a vector of raw residuals \code{(observed - fitted)}.}
#'   \item{link}{ a list with elements "\code{mean}" and "\code{dispersion}"
#'         containing the link objects for the respective submodels.}
#'   \item{vcov}{asymptotic covariance matrix of the maximum likelihood
#'       estimator of the model parameters vector. Specifically, the inverse of
#'       the Fisher information matrix.}
#'   \item{logLik}{log-likelihood of the fitted model.}
#'   \item{freq}{expected frequencies after fitting the BerG regression.}
#'   \item{nobs}{the number of observations in the sample.}
#'   \item{df.null}{ residual degrees of freedom in the null model
#'         (constant mean and dispersion), that is, \code{n - 2}.}
#'   \item{df.residual}{ residual degrees of freedom in the fitted model.}
#'   \item{feasible}{logical; if \code{TRUE}, the estimates satisfy the constraints.}
#'   \item{call}{ the function call.}
#'   \item{formula}{the formula used to specify the model in \code{bergreg}.}
#'   \item{terms}{ a list with elements "mean", "dispersion" and "full"
#'          containing the terms objects for the respective models.}
#'   \item{y}{ the response vector (if \code{y = TRUE}).}
#'   \item{x}{ a list with elements "mean" and "dispersion" containing the
#'       model matrices from the respective models (if \code{X = TRUE}).}
#'  }
#'
#' @details The basic formula is of type \code{y ~ x1 + x2 + ... + xp} which
#'  specifies the model for the mean response only with \code{p} explanatory
#'  variables. Following the syntax of the \code{betareg} package
#'  (Cribari-Neto and Zeileis, 2010), the model for the dispersion index, say in terms of
#'  z1, z2, ..., zk, is specified as \code{y ~ x1 + x2 + ... + xp | z1 + z2 +
#'  ... +zk} using functionalities inherited from package \code{Formula}
#'  (Zeileis and Croissant, 2010).
#'
#' @references Bourguignon, M. & Medeiros, R. (2021). A simple and useful regression model for fitting count data.
#' @references Cribari-Neto F, Zeileis A (2010). Beta regression in R. Journal
#'  of Statistical Software, 34, 1–24
#' @references Zeileis A, Croissant Y (2010). Extended model formulas in
#'  R: multiple parts and multiple responses. Journal of Statistical Software,
#'  34, 1–13.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' # Dataset: grazing
#' data(grazing)
#' head(grazing)
#'
#' layout(matrix(c(1, 2, 1, 3), 2, 2))
#' # Response variable (Number of understorey birds)
#' barplot(table(grazing$birds), xlab = "Number of understorey birds", ylab = "Frequency")
#'
#' # Explanatory variables
#' boxplot(birds ~ when, grazing, ylab = "Number of understorey birds",
#'         xlab = "When the bird count was conduct", pch = 16)
#' boxplot(birds ~ grazed, grazing, ylab = "Number of understorey birds",
#'         xlab = " Which side of the stockproof fence", pch = 16)
#' layout(1)
#'
#' # Fit of the BerG regression model with varying dispersion
#' fit_disp <- bergreg(birds ~ when + grazed | when + grazed, data = grazing)
#' summary(fit_disp)
#'
#' # Fit of the BerG regression model with constant dispersion
#' fit <- bergreg(birds ~ when + grazed, data = grazing, link.phi = "identity")
#' summary(fit)
#'
#' # Diagnostic
#' layout(matrix(c(1, 3, 5, 2, 4, 5), 3, 2))
#' plot(fit, ask = FALSE)
#' plot(fit, which = 3:4, ask = FALSE)
#' plot(fit, which = 5, ask = FALSE)
#' layout(1)
#'
NULL

#' @rdname bergreg
#' @export
bergreg <- function(formula, data, link = c("log", "sqrt", "identity"),
                    link.phi = "log", disp.test = FALSE,
                    optimizer = c("nloptr", "optim"),
                    y = FALSE, x = FALSE,
                    control = berg_control(...), ...)
{
  responseLogical <- y
  xLogical <- x

  ## Call
  cl <- match.call()
  if (missing(data))  data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## Formula
  oformula <- stats::as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if (length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  }else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula

  ## Evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## Extract terms, model matrix, response
  mt <- stats::terms(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2L))
  y <- stats::model.response(mf, "numeric")
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)
  p <- NCOL(X); k <- NCOL(Z); n <- length(y)

  if (length(y) < 1)
    stop("empty model")
  if (any(y < 0))
    stop("invalid dependent variable, all observations must be non-negatives integers")

  if (is.character(link)){
    link <- match.arg(link, c("log", "sqrt", "identity"))
  }

  if (is.character(link.phi)){
    link.phi <- match.arg(link.phi, c("log", "sqrt", "identity"))
  }

  opt <- mle_berg(y, X, Z, link = link, link.phi = link.phi, control = control,
                  optimizer = optimizer)

  ## Coefficients
  beta <- opt$est[1:p]
  gamma <- opt$est[1:k + p]

  ## Covariance matrix
  vcov <- K_berg(opt$est, X, Z, link = link, link.phi = link.phi, inverse = TRUE)

  ## Fitted values
  mu <- as.numeric(g(link)$inv(X%*%beta))
  phi <- as.numeric(g(link.phi)$inv(Z%*%gamma))

  ## Expected frequencies
  expect_berg <- function(y){
    x <- sort(unique(y))
    n <- length(y)
    s <- rep(0, length(x))

    if (length(phi) == 1) phi = as.matrix(rep(phi, n))

    for(i in 1:n){
      s <- s + dberg(x, mu[i], phi[i])
    }
    return(s)
  }

  freq <- expect_berg(y)

  ## set up return value
  out <- list(coefficients = list(mean = beta, dispersion = gamma),
              fitted.values = structure(mu, .Names = names(y)),
              phi = phi,
              residuals = y - mu,
              link = list(mean = link, dispersion = link.phi),
              vcov = vcov,
              logLik = opt$logLik,
              freq = freq,
              nobs = n,
              df.null = n - 2,
              df.residual = n - p - k,
              feasible = opt$feasible)

  ## Further model information
  out$call <- cl
  out$formula <- formula
  out$terms <- list(mean = mtX, dispersion = mtZ, full = mt)
  if(responseLogical) out$y <- y
  if(xLogical) out$x <- list(mean = X, dispersion = Z)


  class(out) <- "bergreg"
  out
}

