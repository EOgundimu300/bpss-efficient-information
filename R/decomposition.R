# Local exclusion decomposition (Section 3.2 of the paper). With the excluded
# covariate residualized on the outcome covariates, the residual gradient along
# the exclusion path is g(gamma_V) = g(0) + gamma_V h + r, where h is the
# analytic path derivative at gamma_V = 0. With n^{-1}-normalized norms in the
# information inner product at the evaluation point,
#   K0 = |g(0)|^2,  K1 = |h|^2,  K01 = <g(0), h>,
#   A = K0 + gamma_V^2 K1,  B = A + 2 gamma_V K01.

set_gamma_component <- function(theta, idx_V, gamma_V) { th <- as.numeric(theta); th[idx_V] <- gamma_V; th }

weighted_gradient_inner <- function(G1, G2, prob_w) {
  G1 <- as.matrix(G1); G2 <- as.matrix(G2); prob_w <- as.matrix(prob_w)
  if (!all(dim(G1) == dim(G2)) || !all(dim(G1) == dim(prob_w)))
    stop("G1, G2 and prob_w must have the same dimensions.")
  sum(G1 * G2 / prob_w)
}

# Analytic derivative h of the residual gradient along the exclusion path at
# gamma_V = 0 (expression (S1) of the Supplementary Material).
bpss_residualized_gradient_derivative <- function(theta0, XS, XO, idx_V,
                                                  eps = 1e-12,
                                                  rcond_tol = 1e-12) {
  XS <- as.matrix(XS); XO <- as.matrix(XO)
  dd <- theta_dims(XS, XO)
  if (length(theta0) != dd$ptheta)
    stop("theta0 length is inconsistent with XS and XO.")
  if (idx_V < 1L || idx_V > dd$nS)
    stop("idx_V must index a selection-equation coefficient.")
  if (abs(theta0[idx_V]) > 100 * .Machine$double.eps)
    stop("theta0 must have gamma_V = 0.")

  info0 <- bpss_info_array(theta0, XS, XO, eps = eps,
                           rcond_tol = rcond_tol)
  if (isTRUE(info0$probability_floor_hit))
    stop("analytic decomposition requires all baseline cell probabilities ",
         "to exceed eps.")

  cd <- info0$cell
  n <- nrow(XS)
  ptheta <- dd$ptheta
  idx_psi <- dd$idx_psi
  idx_eta <- dd$idx_eta
  v <- XS[, idx_V]
  a <- cd$a; b <- cd$b; rho <- cd$rho
  phi_b <- stats::dnorm(b)
  B <- cd$B
  phi2 <- cd$phi2

  ## Second derivatives required for differentiation with respect to
  ## gamma_V. Since db/dgamma_V = v, these are multiplied by v below.
  p0_bb  <- b * phi_b
  p11_bb <- -b * B - rho * phi2
  p10_bb <- -b * phi_b - p11_bb
  p11_ab <- phi2
  p10_ab <- -phi2
  p11_b_eta <- -(b - rho * a) * phi2
  p10_b_eta <- -p11_b_eta

  dotD <- array(0, dim = dim(cd$D), dimnames = dimnames(cd$D))
  idx_gamma <- dd$idx_gamma
  idx_beta <- dd$idx_beta

  dotD[, "p0", idx_gamma] <-
    sweep(XS, 1, p0_bb * v, "*")
  dotD[, "p10", idx_gamma] <-
    sweep(XS, 1, p10_bb * v, "*")
  dotD[, "p11", idx_gamma] <-
    sweep(XS, 1, p11_bb * v, "*")
  dotD[, "p10", idx_beta] <-
    sweep(XO, 1, p10_ab * v, "*")
  dotD[, "p11", idx_beta] <-
    sweep(XO, 1, p11_ab * v, "*")
  dotD[, "p10", idx_eta] <- p10_b_eta * v
  dotD[, "p11", idx_eta] <- p11_b_eta * v

  dot_prob <- cd$D[, , idx_V, drop = TRUE]
  Dflat <- matrix(cd$D, nrow = n * 3L, ncol = ptheta)
  dotDflat <- matrix(dotD, nrow = n * 3L, ncol = ptheta)
  pvec <- as.vector(cd$prob)
  dpvec <- as.vector(dot_prob)

  dotI <- crossprod(dotDflat, Dflat / pvec) +
          crossprod(Dflat, dotDflat / pvec) -
          crossprod(Dflat, Dflat * (dpvec / pvec^2))
  dotI <- (dotI + t(dotI)) / 2

  Ipp <- info0$I[idx_psi, idx_psi, drop = FALSE]
  dotIpp <- dotI[idx_psi, idx_psi, drop = FALSE]
  dotIpe <- dotI[idx_psi, idx_eta, drop = FALSE]
  M0 <- info0$M
  R <- chol(Ipp)
  dotM <- backsolve(R, forwardsolve(t(R), dotIpe - dotIpp %*% M0))

  Deta_dot <- dotD[, , idx_eta, drop = TRUE]
  Dpsi_dot <- matrix(dotD[, , idx_psi, drop = FALSE],
                     nrow = n * 3L, ncol = length(idx_psi))
  Dpsi0 <- matrix(cd$D[, , idx_psi, drop = FALSE],
                  nrow = n * 3L, ncol = length(idx_psi))
  H <- Deta_dot -
       matrix(drop(Dpsi_dot %*% M0 + Dpsi0 %*% dotM),
              nrow = n, ncol = 3L)
  colnames(H) <- c("p0", "p10", "p11")

  list(H = H, dotD = dotD, dot_prob = dot_prob,
       dotI = dotI, dotM = dotM, info0 = info0,
       max_abs_row_sum = max(abs(rowSums(H))))
}

# Centred finite-difference version of h, used only to check the analytic
# derivative (validate_fd = TRUE).
.bpss_fd_residualized_gradient_derivative <- function(
    theta0, XS, XO, idx_V, h = 1e-5,
    eps = 1e-12, rcond_tol = 1e-12) {
  if (!is.finite(h) || h <= 0) stop("h must be positive and finite.")
  theta_plus <- set_gamma_component(theta0, idx_V, h)
  theta_minus <- set_gamma_component(theta0, idx_V, -h)
  Gplus <- bpss_info_array(
    theta_plus, XS, XO, eps = eps,
    rcond_tol = rcond_tol)$residualized_gradient
  Gminus <- bpss_info_array(
    theta_minus, XS, XO, eps = eps,
    rcond_tol = rcond_tol)$residualized_gradient
  (Gplus - Gminus) / (2 * h)
}

# Decomposition at a given internal parameter vector theta; idx_V indexes the
# excluded covariate among the columns of XS. Weights are evaluated at theta;
# the baseline gradient and the path derivative at gamma_V = 0.
bpss_decomp_theta <- function(
    theta, XS, XO, idx_V, eps = 1e-12, rcond_tol = 1e-12,
    residualize = TRUE, K0_tol = 1e-8,
    validate_fd = FALSE, fd_h = c(1e-3, 1e-4, 1e-5),
    fd_tol = 1e-6, fd_stop = TRUE) {
  XS <- as.matrix(XS); XO <- as.matrix(XO)
  dd <- theta_dims(XS, XO)
  if (length(theta) != dd$ptheta)
    stop("theta length inconsistent with XS and XO.")
  if (idx_V < 1L || idx_V > dd$nS)
    stop("idx_V must index a selection-equation coefficient.")
  if (!is.finite(K0_tol) || K0_tol < 0)
    stop("K0_tol must be nonnegative and finite.")
  n <- nrow(XS)
  gamma_V <- unname(theta[idx_V])

  ## Work in the OLS-residualised coordinate system used in the paper.
  alpha <- NULL
  if (isTRUE(residualize)) {
    if (is.null(colnames(XO)) || is.null(colnames(XS)))
      stop("residualization requires column names on XO and XS.")
    if (colnames(XS)[idx_V] %in% colnames(XO))
      stop("the named covariate is not excluded from the outcome equation.")
    pos <- match(colnames(XO), colnames(XS))
    if (anyNA(pos))
      stop("every outcome covariate must appear in the selection equation.")
    v <- XS[, idx_V]
    qx <- qr(XO)
    if (qx$rank < ncol(XO))
      stop("the outcome design is rank deficient.")
    alpha <- as.numeric(qr.coef(qx, v))
    if (anyNA(alpha))
      stop("OLS residualisation of the excluded covariate failed.")
    XS_use <- XS
    XS_use[, idx_V] <- v - as.numeric(XO %*% alpha)
    theta_use <- as.numeric(theta)
    theta_use[pos] <- theta_use[pos] + gamma_V * alpha
  } else {
    XS_use <- XS
    theta_use <- as.numeric(theta)
  }
  names(theta_use) <- make_theta_names(XS_use, XO)

  info_eval <- bpss_info_array(theta_use, XS_use, XO,
                               eps = eps,
                               rcond_tol = rcond_tol)
  theta0 <- set_gamma_component(theta_use, idx_V, 0)
  deriv <- bpss_residualized_gradient_derivative(
    theta0, XS_use, XO, idx_V = idx_V,
    eps = eps, rcond_tol = rcond_tol
  )
  info0 <- deriv$info0
  G0 <- info0$residualized_gradient
  H <- deriv$H
  Wprob <- info_eval$cell$prob_w

  fd_validation <- NULL
  if (isTRUE(validate_fd)) {
    fd_h <- sort(unique(as.numeric(fd_h)))
    fd_h <- fd_h[is.finite(fd_h) & fd_h > 0]
    if (!length(fd_h))
      stop("fd_h must contain at least one positive finite step.")
    fd_errors <- vapply(fd_h, function(hh) {
      Hfd <- .bpss_fd_residualized_gradient_derivative(
        theta0, XS_use, XO, idx_V = idx_V, h = hh,
        eps = eps, rcond_tol = rcond_tol)
      max(abs(Hfd - H))
    }, numeric(1))
    names(fd_errors) <- format(fd_h, scientific = TRUE)
    best_error <- min(fd_errors)
    fd_validation <- list(
      steps = fd_h, max_abs_errors = fd_errors,
      best_max_abs_error = best_error,
      tolerance = fd_tol, pass = is.finite(best_error) &&
        best_error <= fd_tol
    )
    if (isTRUE(fd_stop) && !isTRUE(fd_validation$pass))
      stop("analytic derivative failed centred finite-difference validation: ",
           "best maximum absolute discrepancy = ",
           format(best_error, digits = 6), call. = FALSE)
  }

  K0  <- weighted_gradient_inner(G0, G0, Wprob) / n
  K01 <- weighted_gradient_inner(G0, H, Wprob) / n
  K1  <- weighted_gradient_inner(H, H, Wprob) / n
  denom <- sqrt(pmax(K0 * K1, 0))
  K0_resolved <- is.finite(K0) && K0 > K0_tol
  # cosine of the angle between g(0) and h; kappa carries the sign of gamma_V,
  # so that kappa > 0 when the cross term 2 gamma_V K01 is positive.
  kappa_unsigned <- if (K0_resolved && is.finite(denom) && denom > eps)
    max(-1, min(1, K01 / denom)) else NA_real_
  kappa <- if (K0_resolved && abs(gamma_V) > eps)
    sign(gamma_V) * kappa_unsigned else NA_real_
  omega <- if (K0_resolved && is.finite(K1) && K1 >= 0)
    abs(gamma_V) * sqrt(K1 / K0) else NA_real_
  undefined_reason <- if (!K0_resolved)
    sprintf("K0 = %.6g is not above K0_tol = %.6g", K0, K0_tol)
  else if (!is.finite(denom) || denom <= eps)
    "the cross-term normalization is zero or non-finite"
  else if (abs(gamma_V) <= eps)
    "gamma_V is zero, so kappa is not signed"
  else NA_character_

  A <- K0 + gamma_V^2 * K1
  B <- K0 + 2 * gamma_V * K01 +
       gamma_V^2 * K1
  J_over_n <- info_eval$J / n
  B_approx_error <- abs(J_over_n - B)
  expansion_relative_gap <- B_approx_error /
    max(abs(J_over_n), abs(B), .Machine$double.eps)

  r_resid <- info_eval$residualized_gradient -
             G0 - gamma_V * H
  remainder_norm2 <-
    weighted_gradient_inner(r_resid, r_resid, Wprob) / n

  list(
    idx_V = idx_V, gamma_V = gamma_V, n = n,
    I_etaeta = info_eval$I_etaeta, J = info_eval$J,
    J_over_n = J_over_n,
    J_over_I_etaeta = info_eval$J / info_eval$I_etaeta,
    Delta = info_eval$Delta,
    K0 = K0, K1 = K1, K01 = K01,
    kappa = kappa, kappa_unsigned = kappa_unsigned, omega = omega,
    K0_tol = K0_tol, K0_resolved = K0_resolved,
    kappa_omega_undefined_reason = undefined_reason,
    cross_term = 2 * gamma_V * K01,
    exclusion_term = gamma_V^2 * K1,
    A = A, B = B,
    B_approx_error = B_approx_error,
    expansion_relative_gap = expansion_relative_gap,
    remainder_norm2 = remainder_norm2,
    residualize = isTRUE(residualize), alpha = alpha,
    analytic_derivative = H,
    derivative_row_sum_error = deriv$max_abs_row_sum,
    finite_difference_validation = fd_validation,
    fd_validation_max_diff = if (is.null(fd_validation)) NA_real_ else
      fd_validation$best_max_abs_error,
    fd_validation_pass = if (is.null(fd_validation)) NA else
      fd_validation$pass,
    baseline_information = info0,
    evaluated_information = info_eval
  )
}

#' Local exclusion decomposition of the efficient information
#'
#' Evaluates the components of the local expansion in Section 3.2 of the paper
#' for one excluded covariate \eqn{V}. In the least-squares residualized
#' coordinates, the residual gradient along the exclusion path is expanded at
#' \eqn{\gamma_V = 0}. With norms in the information inner product at the
#' evaluation point and normalized by \eqn{n}, \eqn{K_{0,n}} is the squared norm
#' of the baseline gradient, \eqn{K_{1,n}} that of the path derivative and
#' \eqn{K_{01,n}} their inner product. Then \eqn{A_n = K_{0,n} + \gamma_V^2
#' K_{1,n}}, \eqn{B_n = A_n + 2 \gamma_V K_{01,n}}, and \eqn{\kappa_n} and
#' \eqn{\omega_n} are as defined in Section S4.2. The expansion is local in
#' \eqn{\gamma_V}: \eqn{B_n} need not approximate \eqn{J_n/n} when
#' \eqn{\gamma_V} is not small.
#'
#' \code{bpss_decompose()} evaluates the decomposition at the unrestricted
#' estimate (\code{at = "mle"}) or at the constrained fit under \eqn{\rho = 0}
#' (\code{at = "null"}) of a fitted model. \code{bpss_decomposition()}
#' evaluates it at supplied coefficients and design matrices, as in the
#' design-level calculations of the paper.
#'
#' @param object A fitted object of class \code{"bpss"}.
#' @param exclusion Name of the excluded selection covariate (a column of the
#'   selection design matrix that is not in the outcome equation). For
#'   \code{bpss_decompose()} it may be omitted when exactly one covariate is
#'   excluded; for \code{bpss_decomposition()} it may also be a column index of
#'   \code{XS}.
#' @param at Where to evaluate the decomposition for a fitted model.
#' @param residualize Logical; work in the residualized coordinates of the
#'   paper. Requires every outcome column to appear in the selection design.
#' @param K0_tol \eqn{\kappa_n} and \eqn{\omega_n} are not reported when
#'   \eqn{K_{0,n}} does not exceed this value.
#' @param beta,gamma,rho,XO,XS As in \code{\link{bpss_information}}; the
#'   matrices need column names.
#' @param eps,rcond_tol As in \code{\link{bpss_information}}.
#'
#' @return An object of class \code{"bpss_decomp"}: a list with components
#'   \code{K0}, \code{K1}, \code{K01}, \code{cross_term}
#'   (\eqn{2\gamma_V K_{01,n}}), \code{exclusion_term}
#'   (\eqn{\gamma_V^2 K_{1,n}}), \code{A}, \code{B}, \code{kappa},
#'   \code{omega}, \code{J_over_n}, \code{J_over_I_etaeta}, \code{Delta},
#'   \code{gamma_V}, \code{remainder_norm2} (squared norm of the remainder
#'   \eqn{r_{i,n}}) and further details.
#' @seealso \code{\link{bpss_information}}
#' @export
bpss_decompose <- function(object, exclusion = NULL, at = c("mle", "null"),
                           residualize = TRUE, K0_tol = 1e-8) {
  at <- match.arg(at)
  if (!inherits(object, "bpss"))
    stop("object must be of class 'bpss'.")
  XS <- object$XS; XO <- object$XO; D <- object$D
  if (is.null(exclusion)) {
    if (length(object$excluded) == 1L)
      exclusion <- object$excluded
    else
      stop("name one excluded covariate in 'exclusion'.")
  }
  if (length(exclusion) != 1L)
    stop("'exclusion' must name a single covariate.")
  idx_V <- match(exclusion, colnames(XS))
  if (is.na(idx_V))
    stop("exclusion '", exclusion,
         "' is not a selection-equation covariate.")

  if (at == "mle") {
    if (!isTRUE(object$unrestricted$ok))
      stop("unrestricted fit unavailable; cannot decompose at the estimate.")
    theta <- object$unrestricted$theta
  } else {
    fc <- bpss_fit_null_factorized(
      D, XS, XO, eps = object$control$eps,
      rcond_tol = object$control$rcond_tol,
      coef_bound = object$control$coef_bound,
      maxit = object$control$maxit,
      optim_factr = object$control$optim_factr,
      optim_pgtol = object$control$optim_pgtol,
      stationarity_tol = object$control$stationarity_tol,
      boundary_fraction = object$control$boundary_fraction,
      check_generic = TRUE
    )
    if (!isTRUE(fc$ok))
      stop("constrained fit at rho = 0 failed.")
    theta <- fc$theta
  }
  dec <- bpss_decomp_theta(
    theta, XS, XO, idx_V = idx_V,
    eps = object$control$eps,
    rcond_tol = object$control$rcond_tol,
    residualize = residualize, K0_tol = K0_tol
  )
  dec$exclusion <- exclusion
  dec$at <- at
  class(dec) <- "bpss_decomp"
  dec
}

#' @rdname bpss_decompose
#' @examples
#' n <- 500
#' x <- seq(-1.5, 1.5, length.out = n)
#' XO <- cbind("(Intercept)" = 1, x = x)
#' XS <- cbind("(Intercept)" = 1, x = x, v = x^2 - mean(x^2))
#' bpss_decomposition(beta = c(-0.5, 0.4), gamma = c(0.4, 0.3, 0.2),
#'                    rho = 0.3, XO = XO, XS = XS, exclusion = "v")
#' @export
bpss_decomposition <- function(beta, gamma, rho, XO, XS, exclusion,
                               residualize = TRUE, K0_tol = 1e-8,
                               eps = 1e-12, rcond_tol = 1e-12) {
  XO <- as.matrix(XO); XS <- as.matrix(XS)
  theta <- .bpss_theta_from_coef(beta, gamma, rho, XO, XS)
  idx_V <- if (is.character(exclusion)) match(exclusion, colnames(XS)) else
    as.integer(exclusion)
  if (length(idx_V) != 1L || is.na(idx_V) || idx_V < 1L || idx_V > ncol(XS))
    stop("'exclusion' must name or index one column of XS.")
  dec <- bpss_decomp_theta(theta, XS, XO, idx_V = idx_V, eps = eps,
                           rcond_tol = rcond_tol, residualize = residualize,
                           K0_tol = K0_tol)
  dec$exclusion <- if (!is.null(colnames(XS))) colnames(XS)[idx_V] else
    as.character(idx_V)
  dec$at <- "supplied"
  class(dec) <- "bpss_decomp"
  dec
}
