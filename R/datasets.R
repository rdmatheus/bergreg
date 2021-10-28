#' @title Bird abundance in grazing areas
#'
#' @description The density of understorey birds at a
#' series of sites in two areas either side of a stockproof
#' fence.
#'
#' @format A data frame with 62 observations on the following 3 variables.
#'
#' \itemize{
#' \item \code{birds}: the number of understorey birds; a numeric vector.
#' \item \code{when}: when the bird count was conducted; a factor with levels \code{Before} (before herbivores were removed) and \code{After} (after herbivores were removed).
#' \item \code{grazed}: which side of the stockproof fence; a factor with levels \code{Reference} (grazed by native herbivores) and \code{Feral} (grazed by feral herbivores, mainly horses).
#' }
#'
#' @usage data(grazing)
#'
#' @details In this experiment, the density of understorey birds at
#' a series of sites in two areas either side of a stockproof fence
#' were compared. Once side had limited grazing (mainly from native
#' herbivores), and the other was heavily grazed by feral herbivores,
#' mostly horses. Bird counts were done at the sites either side of
#' the fence (the Before measurements). Then the herbivores were
#' removed, and bird counts done again (the After measurements).
#' The measurements are the total number of understorey-foraging
#' birds observed in three 20-minute surveys of two hectare quadrats.
#'
#' @references Alison L. Howes, Martine Maron and Clive A. McAlpine (2010) Bayesian networks and adaptive management of wildlife habitat. Conservation Biology. 24(4), 974-983.2.
#'
"grazing"


#' @title Takeover bids data
#'
#' @description Data of the number of bids received by 126 U.S.
#' firms that were targets of tender offers from 1978 to 1985, along with some
#' explanatory variables on the defensive actions taken by management of target
#' firm, firm-specific characteristics and intervention taken by federal
#' regulators.
#'
#' @format A data frame with 126 observations on the following 12 variables.
#'
#' \itemize{
#' \item \code{numbids}: number of bids recevied after the initial bid.
#' \item \code{leglrest}: indicator variable for legal defence by lawsuit.
#' \item \code{realrest}: indicator variable for proposed changes in asset structure.
#' \item \code{finrest}: indicator variable for proposed change in ownership structure.
#' \item \code{whtknght}: indicator varible for management invitation for friendly third-party bid.
#' \item \code{bidprem}: bid price divided by price 14 working days before bid.
#' \item \code{insthold}: percentage of stock held by institutions.
#' \item \code{size}: total book value of assets in billions of dollars.
#' \item \code{regulatn}: indicator variable for Department of Justice intervention.
#' }
#'
#' @usage data(bids)
#'
#' @details Data originally from Jaggia and Thosar (1993) and used
#' as an example in Cameron and Trivedi (2013) and Saez-Castillo
#' and Conde-Sanchez (2013).
#'
#' @references Cameron, A.C., Trivedi, P.K. (2013). Regression Analysis of Count Data. Cambridge University Press, second edition. \cr
#'
#' Jaggia, S., Thosar, S. (1993). Multiple Bids as a Consequence of Target Management Resistance. Review of Quantitative Finance and Accounting, 447-457. \cr
#'
#' Saez-Castillo, A.J., Conde-Sanchez, A. (2013). A hyper-Poisson regression model for overdispersed and underdispersed count data. Computational Statistics and Data Analysis, 61, 148-157.doi: 10.1016/j.csda.2012.12.009
#'
"bids"
