#' bpssInfo: efficient information and score inference for the dependence
#' parameter in bivariate probit sample selection models
#'
#' Implements the methods of Ogundimu, "Efficient information for the
#' dependence parameter in bivariate probit sample selection models": the
#' efficient (nuisance-adjusted) information for
#' \eqn{\eta = \mathrm{atanh}(\rho)}{eta = atanh(rho)} and its local
#' decomposition under an exclusion restriction, and inference on
#' \eqn{\rho}{rho} by inverting the score test with the regression coefficients
#' estimated at the hypothesized correlation, together with the relative
#' profile likelihood and the information diagnostics recommended in the paper.
#'
#' The main function is \code{\link{bpss}}. Design-level calculations at
#' specified parameter values use \code{\link{bpss_information}} and
#' \code{\link{bpss_decomposition}}.
#'
#' @keywords internal
#' @importFrom stats coef vcov logLik nobs confint predict fitted
#' @importFrom utils head tail
"_PACKAGE"
