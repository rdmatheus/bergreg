% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/berg_methods.R
\name{predict.bergreg}
\alias{predict.bergreg}
\title{Predict Method for BerG Fits}
\usage{
\method{predict}{bergreg}(
  object,
  newdata = NULL,
  type = c("response", "link", "dispersion", "variance", "quantile"),
  at = 0.5,
  na.action = stats::na.pass,
  ...
)
}
\arguments{
\item{object}{an \code{'bergreg'} object.}

\item{newdata}{optionally, a data frame in which to look for variables
with which to predict. If omitted, the fitted linear predictors are
used.}

\item{type}{the type of prediction required. The default is on the scale of
the response variable \code{("response")}, that is, the fitted values
(fitted means). The alternative \code{"link"} is on the scale of
the linear predictors, that is, predictions are of log-means. The
specification \code{"dispersion"} provides the fitted dispersion
indixes, while \code{"variance"} provides the fitted variances. Finally,
the option \code{"quantile"} gives the fitted quantiles in the order
specified via \code{'at'}.}

\item{at}{the order of the quantile to be predicted if
\code{type = "quantile"}. The default is to predict the median,
that is, \code{at = 0.5}.}

\item{na.action}{function determining what should be done with missing
values in newdata. The default is to predict \code{NA}.}

\item{...}{arguments passed to or from other methods.}
}
\value{
A vector of predictions.
}
\description{
Obtains predictions from a fitted BerG regression object.
}
