# Model primitives for the bivariate probit sample selection model: parameter
# indexing, observed cells, cell probabilities, the observed-data likelihood and
# score, the expected information, its Schur complement J for eta = atanh(rho),
# and the efficient score.
#
# Internal parameter order: theta = (gamma [selection], beta [outcome], eta).
# The S3 methods display (beta, gamma, eta), the order used in the paper.

`%||%` <- function(a, b) if (!is.null(a)) a else b

clamp <- function(x, lower, upper) pmin(pmax(x, lower), upper)

theta_dims <- function(XS, XO) {
  XS <- as.matrix(XS); XO <- as.matrix(XO)
  nS <- ncol(XS); nO <- ncol(XO)
  list(nS = nS, nO = nO, ppsi = nS + nO, ptheta = nS + nO + 1L,
       idx_gamma = seq_len(nS), idx_beta = nS + seq_len(nO),
       idx_eta = nS + nO + 1L, idx_psi = seq_len(nS + nO))
}

make_theta_names <- function(XS, XO) {
  XS <- as.matrix(XS); XO <- as.matrix(XO)
  cnS <- colnames(XS) %||% paste0("s", seq_len(ncol(XS)))
  cnO <- colnames(XO) %||% paste0("o", seq_len(ncol(XO)))
  c(paste0("gamma:", cnS), paste0("beta:", cnO), "eta")
}

make_theta_from_psi_eta <- function(psi, eta, XS, XO) {
  dd <- theta_dims(XS, XO)
  if (length(psi) != dd$ppsi) stop("psi has wrong length.")
  theta <- c(as.numeric(psi), eta = as.numeric(eta))
  names(theta) <- make_theta_names(XS, XO)
  theta
}

# Observed cells: 0 (S = 0), 10 (S = 1, Y = 0), 11 (S = 1, Y = 1).
observed_cell_matrix <- function(S, Y_obs) {
  if (length(S) != length(Y_obs)) stop("S and Y_obs must have the same length.")
  if (anyNA(Y_obs[S == 1L])) stop("Y_obs contains NA among selected observations.")
  Y_fit <- ifelse(S == 1L, Y_obs, 0L)
  out <- cbind(p0  = as.integer(S == 0L),
               p10 = as.integer(S == 1L & Y_fit == 0L),
               p11 = as.integer(S == 1L & Y_fit == 1L))
  storage.mode(out) <- "numeric"
  out
}

D_to_YS_YO <- function(D) {
  if (!all(c("p0", "p10", "p11") %in% colnames(D)))
    stop("D must have columns 'p0', 'p10' and 'p11'.")
  list(YS = as.numeric(D[, "p10"] + D[, "p11"]),
       YO = as.numeric(D[, "p11"]))
}

# Starting values from separate probit fits.
safe_glm_probit_coef <- function(x, y, coef_bound = 8) {
  x <- as.matrix(x); y <- as.numeric(y); p <- ncol(x)
  if (length(unique(y[is.finite(y)])) < 2L) return(rep(0, p))
  fit <- tryCatch(
    suppressWarnings(stats::glm.fit(x = x, y = y,
                                    family = stats::binomial(link = "probit"),
                                    control = stats::glm.control(maxit = 50))),
    error = function(e) NULL)
  if (is.null(fit) || is.null(fit$coefficients)) return(rep(0, p))
  b <- as.numeric(fit$coefficients); b[!is.finite(b)] <- 0
  clamp(b, -coef_bound, coef_bound)
}

make_marginal_start <- function(D, XS, XO, eta_start = 0, coef_bound = 8) {
  XS <- as.matrix(XS); XO <- as.matrix(XO)
  dy <- D_to_YS_YO(D); S <- dy$YS; Y <- dy$YO
  gamma_start <- safe_glm_probit_coef(XS, S, coef_bound = coef_bound)
  if (sum(S == 1L) >= ncol(XO) + 2L && length(unique(Y[S == 1L])) >= 2L) {
    beta_start <- safe_glm_probit_coef(XO[S == 1L, , drop = FALSE], Y[S == 1L], coef_bound = coef_bound)
  } else {
    beta_start <- rep(0, ncol(XO))
  }
  theta <- c(gamma_start, beta_start, eta = eta_start)
  names(theta) <- make_theta_names(XS, XO)
  theta
}

.bpss_phi2_vec <- function(b, a, rho) {
  q <- 1 - rho^2
  if (!is.finite(q) || q <= 0) stop("rho outside (-1,1).")
  exp(-(b^2 - 2 * rho * b * a + a^2) / (2 * q)) /
    (2 * pi * sqrt(q))
}

.bpss_prob_core <- function(a, b, rho) {
  if (length(a) != length(b)) stop("a and b have different lengths.")
  if (!is.finite(rho) || abs(rho) >= 1) stop("rho outside (-1,1).")

  p0_direct <- stats::pnorm(-b)
  ps_direct <- stats::pnorm(b)
  use_p0 <- p0_direct <= ps_direct
  p0 <- ifelse(use_p0, p0_direct, 1 - ps_direct)
  ps <- 1 - p0

  p10_direct <- as.numeric(pbivnorm::pbivnorm(b, -a, -rho))
  p11_direct <- as.numeric(pbivnorm::pbivnorm(b,  a,  rho))
  if (any(!is.finite(p10_direct)) || any(!is.finite(p11_direct)))
    stop("non-finite selected-cell probability.")

  tol <- 100 * .Machine$double.eps
  bad <- p10_direct < -tol | p11_direct < -tol |
    p10_direct > ps + tol | p11_direct > ps + tol
  if (any(bad))
    stop("selected-cell probability outside its admissible range.")

  p10_direct <- pmin(pmax(p10_direct, 0), ps)
  p11_direct <- pmin(pmax(p11_direct, 0), ps)
  use_p10 <- p10_direct <= p11_direct

  p10 <- numeric(length(ps))
  p11 <- numeric(length(ps))
  p10[use_p10] <- p10_direct[use_p10]
  p11[use_p10] <- ps[use_p10] - p10[use_p10]
  p11[!use_p10] <- p11_direct[!use_p10]
  p10[!use_p10] <- ps[!use_p10] - p11[!use_p10]

  p10[p10 < 0 & p10 >= -tol] <- 0
  p11[p11 < 0 & p11 >= -tol] <- 0

  prob <- cbind(p0 = p0, p10 = p10, p11 = p11)
  if (any(!is.finite(prob)) || any(prob < 0))
    stop("non-finite or negative cell probability.")

  list(prob = prob, p0 = p0, p10 = p10, p11 = p11, ps = ps,
       route_error = max(abs(p10 + p11 - ps)))
}

# Cell probabilities (p0, p10, p11) at theta.
bpss_cell_prob <- function(theta, XS, XO) {
  XS <- as.matrix(XS); XO <- as.matrix(XO)
  nS <- ncol(XS); nO <- ncol(XO)
  if (length(theta) != nS + nO + 1L)
    stop("theta length is inconsistent with XS and XO.")
  gamma <- theta[seq_len(nS)]
  beta <- theta[nS + seq_len(nO)]
  eta <- theta[nS + nO + 1L]
  rho <- tanh(eta)
  a <- drop(XO %*% beta)
  b <- drop(XS %*% gamma)
  .bpss_prob_core(a, b, rho)$prob
}

bpss_loglik <- function(theta, XS, XO, D = NULL, S = NULL,
                         Y_obs = NULL, eps = 1e-12) {
  if (is.null(D)) {
    if (is.null(S) || is.null(Y_obs))
      stop("Provide either D or both S and Y_obs.")
    D <- observed_cell_matrix(S, Y_obs)
  }
  D <- as.matrix(D)
  prob <- tryCatch(
    bpss_cell_prob(theta, XS, XO),
    error = function(e) NULL
  )
  if (is.null(prob)) return(-Inf)
  if (!all(dim(D) == dim(prob))) stop("D has incompatible dimensions.")
  obs_prob <- rowSums(D * prob)
  if (any(!is.finite(obs_prob)) || any(obs_prob <= 0)) return(-Inf)
  sum(log(obs_prob))
}

# Observed-data objective (negative average log-likelihood) and its analytic
# gradient. The value is cached, so the paired fn/gr requests of optim() share
# one evaluation. Each row of D must be the indicator vector of the observed
# cell; expected (probability-weighted) log-likelihoods are not supported.
.bpss_make_observed_evaluator <- function(D, XS, XO, eta_fixed = NULL) {
  D <- as.matrix(D); XS <- as.matrix(XS); XO <- as.matrix(XO)
  if (nrow(D) != nrow(XS) || nrow(D) != nrow(XO) || ncol(D) != 3L)
    stop("D, XS and XO have incompatible dimensions.")
  if (any(rowSums(D) != 1)) stop("each row of D must select one observed cell.")

  n <- nrow(D); nS <- ncol(XS); nO <- ncol(XO)
  cell <- max.col(D, ties.method = "first")
  i0 <- cell == 1L; i10 <- cell == 2L; i11 <- cell == 3L
  cache <- new.env(parent = emptyenv())
  cache$par <- NULL
  cache$value <- NULL

  calculate <- function(par) {
    par <- as.numeric(par)
    if (is.null(eta_fixed)) {
      if (length(par) != nS + nO + 1L)
        stop("unrestricted parameter has wrong length.")
      gamma <- par[seq_len(nS)]
      beta <- par[nS + seq_len(nO)]
      eta <- par[nS + nO + 1L]
    } else {
      if (length(par) != nS + nO)
        stop("constrained nuisance parameter has wrong length.")
      gamma <- par[seq_len(nS)]
      beta <- par[nS + seq_len(nO)]
      eta <- eta_fixed
    }
    rho <- tanh(eta)
    a <- drop(XO %*% beta)
    b <- drop(XS %*% gamma)
    pc <- .bpss_prob_core(a, b, rho)

    obs_prob <- numeric(n)
    obs_prob[i0] <- pc$p0[i0]
    obs_prob[i10] <- pc$p10[i10]
    obs_prob[i11] <- pc$p11[i11]
    if (any(!is.finite(obs_prob)) || any(obs_prob <= 0))
      stop("an observed cell has zero or non-finite probability.")

    q <- 1 - rho^2
    sr <- sqrt(pmax(q, .Machine$double.eps))
    phi_a <- stats::dnorm(a)
    phi_b <- stats::dnorm(b)
    A <- phi_a * stats::pnorm((b - rho * a) / sr)
    B <- phi_b * stats::pnorm((a - rho * b) / sr)
    Ceta <- q * .bpss_phi2_vec(b, a, rho)
    D10b <- phi_b - B

    sb <- sa <- se <- numeric(n)
    sb[i0] <- -phi_b[i0] / pc$p0[i0]
    sb[i10] <- D10b[i10] / pc$p10[i10]
    sb[i11] <- B[i11] / pc$p11[i11]
    sa[i10] <- -A[i10] / pc$p10[i10]
    sa[i11] <-  A[i11] / pc$p11[i11]
    se[i10] <- -Ceta[i10] / pc$p10[i10]
    se[i11] <-  Ceta[i11] / pc$p11[i11]

    score_psi <- c(drop(crossprod(XS, sb)), drop(crossprod(XO, sa)))
    score <- if (is.null(eta_fixed)) c(score_psi, sum(se)) else score_psi
    list(
      objective = -sum(log(obs_prob)) / n,
      gradient = -score / n,
      logLik = sum(log(obs_prob)),
      min_observed_probability = min(obs_prob)
    )
  }

  get_value <- function(par) {
    p <- as.numeric(par)
    same <- !is.null(cache$par) && length(cache$par) == length(p) &&
      all(cache$par == p)
    if (!same) {
      cache$par <- p
      cache$value <- calculate(p)
    }
    cache$value
  }

  list(
    fn = function(par) get_value(par)$objective,
    gr = function(par) get_value(par)$gradient,
    value = get_value,
    cache = cache
  )
}

# Cell probabilities and the n x 3 x p array of their derivatives.
bpss_cell_derivatives_vec <- function(theta, XS, XO, eps = 1e-12) {
  XS <- as.matrix(XS); XO <- as.matrix(XO)
  n <- nrow(XS); nS <- ncol(XS); nO <- ncol(XO)
  ptheta <- nS + nO + 1L
  if (nrow(XO) != n) stop("XS and XO must have the same number of rows.")
  if (length(theta) != ptheta) stop("theta length is inconsistent with XS and XO.")
  gamma <- theta[seq_len(nS)]
  beta  <- theta[nS + seq_len(nO)]
  eta   <- theta[ptheta]
  rho <- tanh(eta)
  if (!is.finite(rho) || abs(rho) >= 1) stop("rho outside (-1, 1).")

  a <- drop(XO %*% beta)
  b <- drop(XS %*% gamma)
  prob <- bpss_cell_prob(theta, XS, XO)
  prob_w <- pmax(prob, eps)
  floor_hit <- prob < eps

  s_rho <- sqrt(pmax(1 - rho^2, .Machine$double.eps))
  phi_a <- stats::dnorm(a)
  phi_b <- stats::dnorm(b)
  A <- phi_a * stats::pnorm((b - rho * a) / s_rho)  # d p11 / d a
  B <- phi_b * stats::pnorm((a - rho * b) / s_rho)  # d p11 / d b
  D_b_p10 <- phi_b - B
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  phi2 <- mvtnorm::dmvnorm(cbind(b, a), mean = c(0, 0), sigma = Sigma)
  C_eta <- (1 - rho^2) * phi2                        # d p11 / d eta

  Darr <- array(0, dim = c(n, 3L, ptheta),
                dimnames = list(NULL, c("p0", "p10", "p11"),
                                make_theta_names(XS, XO)))
  idx_gamma <- seq_len(nS)
  idx_beta  <- nS + seq_len(nO)
  idx_eta   <- ptheta

  Darr[, "p0",  idx_gamma] <- sweep(XS, 1, -phi_b, "*")
  Darr[, "p10", idx_gamma] <- sweep(XS, 1, D_b_p10, "*")
  Darr[, "p11", idx_gamma] <- sweep(XS, 1, B, "*")
  Darr[, "p10", idx_beta]  <- sweep(XO, 1, -A, "*")
  Darr[, "p11", idx_beta]  <- sweep(XO, 1, A, "*")
  Darr[, "p10", idx_eta]   <- -C_eta
  Darr[, "p11", idx_eta]   <-  C_eta

  list(prob = prob, prob_w = prob_w, probability_floor_hit = floor_hit,
       D = Darr, rho = rho, eta = eta, a = a, b = b,
       A = A, B = B, phi2 = phi2)
}

# Expected information from the full n x 3 x p derivative array. The exclusion
# decomposition uses this form because its analytic path derivative needs the
# array. Algebraically identical to bpss_info_compact().
bpss_info_array <- function(theta, XS, XO, eps = 1e-12,
                            rcond_tol = 1e-12) {
  cd <- bpss_cell_derivatives_vec(theta, XS, XO, eps = eps)
  Darr <- cd$D
  prob_w <- cd$prob_w
  n <- dim(Darr)[1]
  ptheta <- dim(Darr)[3]
  idx_eta <- ptheta
  idx_psi <- seq_len(ptheta - 1L)

  Dflat <- matrix(Darr, nrow = n * 3L, ncol = ptheta)
  w <- as.vector(prob_w)
  WD <- Dflat / sqrt(w)
  I <- crossprod(WD)
  I <- (I + t(I)) / 2
  dimnames(I) <- list(dimnames(Darr)[[3]], dimnames(Darr)[[3]])

  I_pp <- I[idx_psi, idx_psi, drop = FALSE]
  I_pe <- I[idx_psi, idx_eta, drop = FALSE]
  I_ee <- I[idx_eta, idx_eta]
  rc <- rcond(I_pp)
  if (!is.finite(rc) || rc < rcond_tol)
    stop("I_psipsi is numerically singular.")

  R <- tryCatch(chol(I_pp), error = function(e) NULL)
  if (is.null(R)) stop("I_psipsi is not positive definite.")
  M <- backsolve(R, forwardsolve(t(R), I_pe))

  Deta <- Darr[, , idx_eta, drop = TRUE]
  Dpsi_flat <- matrix(Darr[, , idx_psi, drop = FALSE],
                      nrow = n * 3L, ncol = length(idx_psi))
  proj <- matrix(drop(Dpsi_flat %*% M), nrow = n, ncol = 3L)
  G <- Deta - proj
  colnames(G) <- c("p0", "p10", "p11")

  ## The residual-norm representation is numerically preferable when J is
  ## small; the Schur difference is retained as an identity check.
  v <- rowSums(G^2 / prob_w)
  J_residual <- sum(v)
  J_schur <- as.numeric(I_ee - crossprod(I_pe, M))
  J <- J_residual
  Delta <- if (is.finite(J) && J > 0) max(v) / J else NA_real_

  list(I = I, I_etaeta = I_ee, J = J, J_schur = J_schur,
       J_identity_error = J_residual - J_schur,
       M = M, residualized_gradient = G, v = v, Delta = Delta,
       rcond_psipsi = rc,
       probability_floor_hit = any(cd$probability_floor_hit),
       min_cell_probability = min(cd$prob),
       cell = cd)
}

# Expected information I, its Schur complement J for eta, the residual
# cell-probability gradients g_i, the contributions j_i and Delta_n, computed
# from weighted cross-products of the index derivatives.
bpss_info_compact <- function(theta, XS, XO, eps = 1e-12,
                              rcond_tol = 1e-12) {
  XS <- as.matrix(XS); XO <- as.matrix(XO)
  dd <- theta_dims(XS, XO)
  if (length(theta) != dd$ptheta)
    stop("theta length is inconsistent with XS and XO.")
  gamma <- theta[dd$idx_gamma]
  beta <- theta[dd$idx_beta]
  eta <- theta[dd$idx_eta]
  rho <- tanh(eta)
  a <- drop(XO %*% beta)
  b <- drop(XS %*% gamma)
  pc <- .bpss_prob_core(a, b, rho)
  prob <- pc$prob
  prob_w <- pmax(prob, eps)
  p0 <- prob_w[, "p0"]
  p10 <- prob_w[, "p10"]
  p11 <- prob_w[, "p11"]

  q <- 1 - rho^2
  sr <- sqrt(pmax(q, .Machine$double.eps))
  phi_a <- stats::dnorm(a)
  phi_b <- stats::dnorm(b)
  A <- phi_a * stats::pnorm((b - rho * a) / sr)
  B <- phi_b * stats::pnorm((a - rho * b) / sr)
  phi2 <- .bpss_phi2_vec(b, a, rho)
  Ceta <- q * phi2
  D10b <- phi_b - B

  wgg <- phi_b^2 / p0 + D10b^2 / p10 + B^2 / p11
  wbb <- A^2 * (1 / p10 + 1 / p11)
  common <- B / p11 - D10b / p10
  wgb <- A * common
  wge <- Ceta * common
  wbe <- A * Ceta * (1 / p10 + 1 / p11)
  wee <- Ceta^2 * (1 / p10 + 1 / p11)

  I <- matrix(0, dd$ptheta, dd$ptheta)
  I[dd$idx_gamma, dd$idx_gamma] <- crossprod(XS * wgg, XS)
  I[dd$idx_beta, dd$idx_beta] <- crossprod(XO * wbb, XO)
  I[dd$idx_gamma, dd$idx_beta] <- crossprod(XS * wgb, XO)
  I[dd$idx_beta, dd$idx_gamma] <- t(I[dd$idx_gamma, dd$idx_beta, drop = FALSE])
  I[dd$idx_gamma, dd$idx_eta] <- drop(crossprod(XS, wge))
  I[dd$idx_eta, dd$idx_gamma] <- I[dd$idx_gamma, dd$idx_eta]
  I[dd$idx_beta, dd$idx_eta] <- drop(crossprod(XO, wbe))
  I[dd$idx_eta, dd$idx_beta] <- I[dd$idx_beta, dd$idx_eta]
  I[dd$idx_eta, dd$idx_eta] <- sum(wee)
  I <- (I + t(I)) / 2
  nm <- make_theta_names(XS, XO)
  dimnames(I) <- list(nm, nm)

  Ipp <- I[dd$idx_psi, dd$idx_psi, drop = FALSE]
  Ipe <- I[dd$idx_psi, dd$idx_eta, drop = FALSE]
  Iee <- I[dd$idx_eta, dd$idx_eta]
  rc <- rcond(Ipp)
  if (!is.finite(rc) || rc < rcond_tol)
    stop("I_psipsi is numerically singular.")
  R <- tryCatch(chol(Ipp), error = function(e) NULL)
  if (is.null(R)) stop("I_psipsi is not positive definite.")
  M <- backsolve(R, forwardsolve(t(R), Ipe))

  Mg <- M[seq_len(dd$nS)]
  Mb <- M[dd$nS + seq_len(dd$nO)]
  mg <- drop(XS %*% Mg)
  mb <- drop(XO %*% Mb)
  G <- cbind(
    p0 = phi_b * mg,
    p10 = -Ceta - D10b * mg + A * mb,
    p11 =  Ceta - B * mg - A * mb
  )
  v <- rowSums(G^2 / prob_w)
  J <- sum(v)
  J_schur <- as.numeric(Iee - crossprod(Ipe, M))
  Delta <- if (is.finite(J) && J > 0) max(v) / J else NA_real_

  list(
    I = I, I_etaeta = Iee, J = J, J_schur = J_schur,
    J_identity_error = J - J_schur,
    M = M, residualized_gradient = G, v = v, Delta = Delta,
    rcond_psipsi = rc,
    probability_floor_hit = any(prob < eps),
    min_cell_probability = min(prob),
    cell = list(
      prob = prob, prob_w = prob_w, rho = rho, eta = eta,
      a = a, b = b, A = A, B = B, phi2 = phi2,
      C_eta = Ceta, D_b_p10 = D10b
    )
  )
}

bpss_info_safe <- function(theta, XS, XO, eps = 1e-12,
                           rcond_tol = 1e-12) {
  tryCatch(
    bpss_info_compact(theta, XS, XO, eps = eps,
                                rcond_tol = rcond_tol),
    error = function(e) NULL
  )
}

# Efficient score U, its fitted variance J, the statistic T = U^2 / J and the
# per-observation quantities q_i and j_i (v_i) at theta.
bpss_efficient_score_stat <- function(theta, XS, XO, D,
                                      eps = 1e-12,
                                      rcond_tol = 1e-12) {
  info <- bpss_info_compact(theta, XS, XO, eps = eps,
                                      rcond_tol = rcond_tol)
  prob <- as.matrix(info$cell$prob)
  G <- as.matrix(info$residualized_gradient)
  D <- as.matrix(D)
  if (!all(dim(D) == dim(prob))) stop("D has incompatible dimensions.")
  if (isTRUE(info$probability_floor_hit))
    stop(sprintf("a fitted cell probability is below %g", eps))
  if (any(prob <= 0)) stop("zero fitted cell probability.")

  raw_s <- G / prob
  mean_s <- rowSums(prob * raw_s) / rowSums(prob)
  s_mat <- raw_s - mean_s
  q_i <- rowSums(D * s_mat)
  U <- sum(q_i)
  v_i <- rowSums(prob * s_mat^2)
  J <- sum(v_i)
  T_score <- if (is.finite(J) && J > 0) U^2 / J else NA_real_
  Delta <- if (is.finite(J) && J > 0) max(v_i) / J else NA_real_
  list(
    U = U, J = J, J_info = info$J,
    J_identity_error = J - info$J,
    I_etaeta = info$I_etaeta, T_score = T_score,
    Delta = Delta, q_i = q_i, v_i = v_i, s = s_mat, info = info
  )
}

#' Efficient information for the dependence parameter
#'
#' Computes, at specified coefficients and correlation, the expected
#' information and the efficient information
#' \deqn{J_n = I_{\eta\eta,n} - I_{\eta\psi,n} I_{\psi\psi,n}^{-1} I_{\psi\eta,n}}{J = I_etaeta - I_etapsi I_psipsi^(-1) I_psieta}
#' for \eqn{\eta = \mathrm{atanh}(\rho)}{eta = atanh(rho)}, with \eqn{\psi}{psi}
#' the regression coefficients. Also returned are the retained fraction
#' \eqn{J_n / I_{\eta\eta,n}}{J / I_etaeta} and \eqn{\Delta_n}{Delta}, the largest
#' proportion of \eqn{J_n}{J} contributed by a single observation.
#'
#' @param beta Outcome coefficients, one per column of \code{XO}.
#' @param gamma Selection coefficients, one per column of \code{XS}.
#' @param rho Correlation, a single value in \eqn{(-1, 1)}{(-1, 1)}.
#' @param XO Outcome design matrix, including the intercept column.
#' @param XS Selection design matrix, including the intercept column; its rows
#'   correspond to those of \code{XO}.
#' @param eps Smallest admissible cell probability. The information is not
#'   computed if any cell probability is below \code{eps}. Set \code{eps = 0} to
#'   compute it with the exact probabilities, for example at parameter values
#'   where some cell probabilities are smaller than \eqn{10^{-12}}{1e-12}.
#' @param rcond_tol Smallest admissible reciprocal condition number of the
#'   nuisance information \eqn{I_{\psi\psi,n}}{I_psipsi}.
#'
#' @return A list with components \code{J} (efficient information for
#'   \eqn{\eta}{eta}), \code{J_rho} (\eqn{J_n / (1 - \rho^2)^2}{J / (1 - rho^2)^2},
#'   the efficient information for \eqn{\rho}{rho}), \code{I_etaeta},
#'   \code{J_over_I_etaeta}, \code{Delta}, \code{contributions} (the
#'   per-observation contributions to \eqn{J_n}{J}), \code{I} (the full expected
#'   information, ordered as selection coefficients, outcome coefficients,
#'   \eqn{\eta}{eta}), \code{M} (\eqn{I_{\psi\psi,n}^{-1} I_{\psi\eta,n}}{I_psipsi^(-1) I_psieta}),
#'   \code{residual_gradient} (the \eqn{n \times 3}{n x 3} matrix of residual
#'   cell-probability gradients), \code{min_cell_probability} and
#'   \code{rcond_psipsi}.
#'
#' @examples
#' n <- 500
#' x <- seq(-1.5, 1.5, length.out = n)
#' v <- sin(3 * x)
#' XO <- cbind("(Intercept)" = 1, x = x)
#' XS <- cbind("(Intercept)" = 1, x = x, v = v)
#' info <- bpss_information(beta = c(-0.5, 0.4), gamma = c(0.4, -1, 0.3),
#'                          rho = 0.5, XO = XO, XS = XS)
#' c(J = info$J, ratio = info$J_over_I_etaeta, Delta = info$Delta)
#' @export
bpss_information <- function(beta, gamma, rho, XO, XS, eps = 1e-12,
                             rcond_tol = 1e-12) {
  XO <- as.matrix(XO); XS <- as.matrix(XS)
  theta <- .bpss_theta_from_coef(beta, gamma, rho, XO, XS)
  info <- bpss_info_compact(theta, XS, XO, eps = eps, rcond_tol = rcond_tol)
  if (isTRUE(info$probability_floor_hit))
    stop("a cell probability is below eps = ", format(eps),
         "; set eps = 0 to compute the information with the exact probabilities.")
  list(J = info$J, J_rho = info$J / (1 - rho^2)^2,
       I_etaeta = info$I_etaeta, J_over_I_etaeta = info$J / info$I_etaeta,
       Delta = info$Delta, contributions = info$v, I = info$I, M = info$M,
       residual_gradient = info$residualized_gradient,
       min_cell_probability = info$min_cell_probability,
       rcond_psipsi = info$rcond_psipsi)
}

# Internal parameter vector from user-facing (beta, gamma, rho), with checks.
.bpss_theta_from_coef <- function(beta, gamma, rho, XO, XS) {
  if (nrow(XO) != nrow(XS))
    stop("XO and XS must have the same number of rows.")
  if (length(beta) != ncol(XO) || any(!is.finite(beta)))
    stop("beta must be finite with one value per column of XO.")
  if (length(gamma) != ncol(XS) || any(!is.finite(gamma)))
    stop("gamma must be finite with one value per column of XS.")
  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1)
    stop("rho must be a single value in (-1, 1).")
  theta <- c(as.numeric(gamma), as.numeric(beta), atanh(rho))
  names(theta) <- make_theta_names(XS, XO)
  theta
}
