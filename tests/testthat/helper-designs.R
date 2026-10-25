# Deterministic designs of Sections S5.1 and S5.2 of the Supplementary
# Material, as used in the simulations of the paper.
.std_rms <- function(x, name = "vector") {
  x <- as.numeric(x)
  if (any(!is.finite(x))) stop(name, " contains non-finite values.")
  x <- x - mean(x)
  s <- sqrt(mean(x^2))
  if (!is.finite(s) || s <= .Machine$double.eps) {
    stop(name, " has zero empirical variance.")
  }
  x / s
}

.residualize_on <- function(v, Z, name = "vector") {
  v <- as.numeric(v)
  Z <- as.matrix(Z)
  qz <- qr(Z)
  if (qz$rank < ncol(Z)) {
    stop("Residualisation design is rank deficient for ", name, ".")
  }
  co <- qr.coef(qz, v)
  if (anyNA(co)) stop("Residualisation failed for ", name, ".")
  as.numeric(v - Z %*% co)
}

van_der_corput <- function(n, base = 2L, start = 1L) {
  n <- as.integer(n)
  base <- as.integer(base)
  start <- as.integer(start)
  if (n <= 0L) stop("n must be positive.")
  if (base < 2L) stop("base must be at least 2.")
  k <- start + seq_len(n) - 1L
  out <- numeric(n)
  for (j in seq_along(k)) {
    x <- k[j]
    f <- 1 / base
    r <- 0
    while (x > 0) {
      r <- r + f * (x %% base)
      x <- floor(x / base)
      f <- f / base
    }
    out[j] <- r
  }
  out
}

.safe_uniroot <- function(f, lower = -12, upper = 12, max_expand = 8L,
                          tol = 1e-10, maxiter = 1000L, ...) {
  lo <- lower
  hi <- upper
  flo <- f(lo)
  fhi <- f(hi)
  k <- 0L
  while (is.finite(flo) && is.finite(fhi) && flo * fhi > 0 && k < max_expand) {
    lo <- 1.5 * lo
    hi <- 1.5 * hi
    flo <- f(lo)
    fhi <- f(hi)
    k <- k + 1L
  }
  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) {
    stop(
      "uniroot bracket failed: f(lower)=", signif(flo, 4),
      ", f(upper)=", signif(fhi, 4), "."
    )
  }
  stats::uniroot(
    f, lower = lo, upper = hi, f.lower = flo, f.upper = fhi,
    tol = tol, maxiter = as.integer(maxiter), check.conv = TRUE, ...
  )$root
}

make_design <- function(n, r, beta1, gammax, gammaV, rho,
                        eps_x = 0.05, target_y = 0.30, target_s = 0.65,
                        aux_base = 2L, aux_start = 1L) {
  n <- as.integer(n)
  if (!is.finite(n) || n < 10L) stop("n must be at least 10.")
  if (!is.finite(r) || abs(r) > 1) stop("r must lie in [-1,1].")
  if (!is.finite(rho) || abs(rho) >= 1) stop("rho must lie in (-1,1).")
  if (!is.finite(eps_x) || eps_x <= 0 || eps_x >= 0.5) {
    stop("eps_x must lie in (0,0.5).")
  }
  if (target_y <= 0 || target_y >= 1 || target_s <= 0 || target_s >= 1) {
    stop("target_y and target_s must lie in (0,1).")
  }

  p <- eps_x + (1 - 2 * eps_x) * ((seq_len(n) - 0.5) / n)
  x <- .std_rms(stats::qnorm(p), "x")
  u0 <- van_der_corput(n, base = aux_base, start = aux_start)
  u0 <- eps_x + (1 - 2 * eps_x) * u0
  u0 <- stats::qnorm(u0)
  Z <- cbind("(Intercept)" = 1, x = x)
  u <- .std_rms(.residualize_on(u0, Z, "auxiliary u"), "residualised auxiliary u")
  V <- .std_rms(r * x + sqrt(pmax(0, 1 - r^2)) * u, "excluded covariate V")

  XO <- cbind("(Intercept)" = 1, x = x)
  XS <- cbind("(Intercept)" = 1, x = x, V = V)
  beta0 <- .safe_uniroot(
    function(b0) mean(stats::pnorm(b0 + beta1 * x)) - target_y
  )
  gamma0 <- .safe_uniroot(
    function(g0) mean(stats::pnorm(g0 + gammax * x + gammaV * V)) - target_s
  )
  theta <- c(gamma0, gammax, gammaV, beta0, beta1, atanh(rho))
  names(theta) <- bpssInfo:::make_theta_names(XS, XO)
  list(
    n = n, theta = theta, XS = XS, XO = XO, x = x, V = V,
    corr_emp = stats::cor(x, V),
    params = list(
      r = r, beta0 = beta0, beta1 = beta1,
      gamma0 = gamma0, gammax = gammax, gammaV = gammaV,
      rho = rho, eps_x = eps_x, target_y = target_y,
      target_s = target_s, aux_base = aux_base, aux_start = aux_start
    )
  )
}

paper_designs <- data.frame(
  label = c("reference", "informative_exclusion", "low_information",
            "matched_small_ratio", "matched_large_ratio", "negative_correlation"),
  rho0 = c(0.5, 0.3, 0.5, 0.3, 0.5, -0.5),
  beta1 = c(0.4, 0.4, 0.8, 0.8, 0.4, 0.4),
  gammax = c(-1.0, -0.5, 0.5, 0.5, 0.5, 0.5),
  gammaV = c(0.30, 1.00, 0.15, 0.30, 1.00, 0.15),
  r = c(0.8, 0.4, 0.8, -0.8, -0.4, 0.8),
  stringsAsFactors = FALSE
)
paper_design <- function(k, n) {
  p <- paper_designs[k, ]
  make_design(n, r = p$r, beta1 = p$beta1, gammax = p$gammax,
              gammaV = p$gammaV, rho = p$rho0)
}
