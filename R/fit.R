# Bounded optimisation and the fitters: the unrestricted fit, the constrained
# fit at a fixed eta, their multi-start versions, and the factorised fit at
# eta = 0 (separate selection and selected-outcome probits).

.bpss_projected_gradient <- function(par, gradient, lower, upper,
                                     active_tol = 1e-7) {
  par <- as.numeric(par); gradient <- as.numeric(gradient)
  lower <- as.numeric(lower); upper <- as.numeric(upper)
  if (!all(lengths(list(par, gradient, lower, upper)) == length(par)))
    stop("par, gradient, lower and upper must have equal lengths.")
  scale <- pmax(1, abs(lower), abs(upper), abs(par))
  at_lower <- par <= lower + active_tol * scale
  at_upper <- par >= upper - active_tol * scale
  pg <- gradient
  pg[at_lower & gradient > 0] <- 0
  pg[at_upper & gradient < 0] <- 0
  list(vector = pg,
       norm = if (length(pg)) max(abs(pg)) else 0,
       at_lower = at_lower, at_upper = at_upper)
}

.bpss_boundary_flags <- function(par, lower, upper, eta_index = NULL,
                                 boundary_fraction = 0.99,
                                 active_tol = 1e-7) {
  par <- as.numeric(par); lower <- as.numeric(lower); upper <- as.numeric(upper)
  if (!is.finite(boundary_fraction) || boundary_fraction <= 0 ||
      boundary_fraction >= 1)
    stop("boundary_fraction must lie in (0,1).")
  if (!is.finite(active_tol) || active_tol < 0)
    stop("active_tol must be nonnegative and finite.")
  threshold <- boundary_fraction * pmax(abs(lower), abs(upper))
  hit <- abs(par) >= threshold
  scale <- pmax(1, abs(lower), abs(upper), abs(par))
  at_lower <- par <= lower + active_tol * scale
  at_upper <- par >= upper - active_tol * scale
  at_bound <- at_lower | at_upper
  eta_hit <- eta_at_bound <- FALSE
  if (!is.null(eta_index)) {
    eta_hit <- any(hit[eta_index])
    eta_at_bound <- any(at_bound[eta_index])
  }
  coef_hit <- if (is.null(eta_index)) any(hit) else any(hit[-eta_index])
  coef_at_bound <- if (is.null(eta_index)) any(at_bound) else
    any(at_bound[-eta_index])
  list(hit = hit, boundary_hit = any(hit),
       eta_boundary_hit = eta_hit, coef_boundary_hit = coef_hit,
       at_bound = at_bound, any_at_bound = any(at_bound),
       eta_at_bound = eta_at_bound, coef_at_bound = coef_at_bound,
       at_lower = at_lower, at_upper = at_upper,
       distance_to_lower = par - lower, distance_to_upper = upper - par)
}

.fit_bpss_optim <- function(par0, objective, gradient = NULL,
                             maxit = 1000, coef_bound = 12,
                             eta_bound = atanh(0.995), factr = 1e7,
                             pgtol = 1e-8, stationarity_tol = 1e-4,
                             active_tol = 1e-7, lower = NULL, upper = NULL) {
  par0 <- as.numeric(par0)
  np <- length(par0)
  if (is.null(lower) || is.null(upper)) {
    lower <- rep(-coef_bound, np)
    upper <- rep(coef_bound, np)
    if (np > 1L) {
      lower[np] <- -eta_bound
      upper[np] <- eta_bound
    }
  }
  lower <- as.numeric(lower); upper <- as.numeric(upper)
  if (length(lower) != np || length(upper) != np || any(lower >= upper))
    stop("invalid optimisation bounds.")

  start_raw <- par0
  par0 <- pmin(pmax(par0, lower), upper)
  start_projected <- any(abs(par0 - start_raw) > 0)

  fn <- function(par) {
    if (any(!is.finite(par)) || any(par < lower) || any(par > upper))
      return(.Machine$double.xmax / 100)
    val <- objective(par)
    if (is.finite(val)) val else .Machine$double.xmax / 100
  }
  gr <- if (is.null(gradient)) NULL else function(par) {
    ans <- as.numeric(gradient(par))
    if (length(ans) != np || any(!is.finite(ans)))
      stop("non-finite analytical gradient.")
    ans
  }

  run_one <- function(start, factr_i, pgtol_i, maxit_i) {
    tryCatch(
      stats::optim(start, fn, gr, method = "L-BFGS-B",
                   lower = lower, upper = upper,
                   control = list(maxit = maxit_i, factr = factr_i,
                                  pgtol = pgtol_i)),
      error = function(e) e
    )
  }

  fit1 <- run_one(par0, factr, pgtol, maxit)
  attempts <- list(fit1)
  retry_start <- if (!inherits(fit1, "error") &&
                     all(is.finite(fit1$par))) fit1$par else par0
  need_retry <- inherits(fit1, "error") ||
                !is.finite(fit1$value) || fit1$convergence != 0L
  if (!need_retry && !is.null(gr)) {
    g1 <- tryCatch(gr(fit1$par), error = function(e) rep(NA_real_, np))
    pg1 <- if (all(is.finite(g1)))
      .bpss_projected_gradient(fit1$par, g1, lower, upper,
                               active_tol = active_tol)$norm else Inf
    need_retry <- !is.finite(pg1) || pg1 > stationarity_tol
  }
  if (need_retry) {
    fit2 <- run_one(retry_start, min(factr, 1e4),
                    min(pgtol, stationarity_tol / 10), 2L * maxit)
    attempts[[2L]] <- fit2
  }

  eligible <- vapply(attempts, function(z)
    !inherits(z, "error") && is.finite(z$value) &&
      length(z$par) == np && all(is.finite(z$par)) &&
      all(z$par >= lower) && all(z$par <= upper), logical(1))
  if (!any(eligible)) {
    msg <- paste(vapply(attempts, function(z)
      if (inherits(z, "error")) conditionMessage(z) else
        paste0("convergence code ", z$convergence), character(1)),
      collapse = "; ")
    return(list(ok = FALSE, optim_converged = FALSE, stationary = FALSE,
                feasible = FALSE, error = msg, par = par0,
                value = fn(par0), convergence = NA_integer_,
                projected_gradient_norm = Inf, lower = lower, upper = upper,
                start = start_raw, start_used = par0,
                start_projected = start_projected, attempts = attempts))
  }

  candidates <- attempts[eligible]
  vals <- vapply(candidates, `[[`, numeric(1), "value")
  fit <- candidates[[which.min(vals)]]
  grad_final <- if (is.null(gr)) rep(NA_real_, np) else
    tryCatch(gr(fit$par), error = function(e) rep(NA_real_, np))
  pg <- if (all(is.finite(grad_final)))
    .bpss_projected_gradient(fit$par, grad_final, lower, upper,
                             active_tol = active_tol) else
    list(vector = rep(NA_real_, np), norm = Inf,
         at_lower = rep(FALSE, np), at_upper = rep(FALSE, np))
  optim_converged <- identical(fit$convergence, 0L)
  stationary <- is.finite(pg$norm) && pg$norm <= stationarity_tol
  feasible <- all(fit$par >= lower) && all(fit$par <= upper)

  list(ok = optim_converged && stationary && feasible &&
             is.finite(fit$value),
       optim_converged = optim_converged, stationary = stationary,
       feasible = feasible, par = fit$par, value = fit$value,
       convergence = fit$convergence, message = fit$message %||% "",
       counts = fit$counts, gradient = grad_final,
       projected_gradient = pg$vector,
       projected_gradient_norm = pg$norm,
       at_lower = pg$at_lower, at_upper = pg$at_upper,
       lower = lower, upper = upper,
       start = start_raw, start_used = par0,
       start_projected = start_projected,
       attempts = attempts, optim = fit)
}

make_failed_fit <- function(status, message = NA_character_, theta_length,
                            psi_length = NULL, theta_attempted = NULL,
                            psi_attempted = NULL, eta_star = NA_real_,
                            optim = NULL) {
  list(ok = FALSE, optim_ok = FALSE, likelihood_ok = FALSE,
       info_ok = FALSE, score_ok = FALSE, wald_usable = FALSE,
       score_usable = FALSE, status = status, message = message,
       eta_star = eta_star, theta = rep(NA_real_, theta_length),
       psi = if (is.null(psi_length)) NULL else rep(NA_real_, psi_length),
       theta_attempted = theta_attempted, psi_attempted = psi_attempted,
       logLik = NA_real_, eta_hat = NA_real_, rho_hat = NA_real_,
       J = NA_real_, I_etaeta = NA_real_, Delta = NA_real_, se_eta = NA_real_,
       U = NA_real_, T_score = NA_real_,
       boundary_hit = NA, eta_boundary_hit = NA,
       coef_boundary_hit = NA, any_at_bound = NA,
       eta_at_bound = NA, coef_at_bound = NA,
       projected_gradient_norm = NA_real_,
       projected_gradient_norm_raw = NA_real_,
       stationary = FALSE, optim_converged = FALSE, feasible = FALSE,
       audit_ok = FALSE, reference_dominated = NA,
       default_discordance = NA, default_logLik_gain = NA_real_,
       default_parameter_gap = NA_real_,
       multiple_near_maxima = NA, branch_ambiguous = NA,
       start_discordance = NA, optim = optim)
}

.bpss_unique_starts <- function(starts, tol = 1e-10) {
  starts <- lapply(starts, as.numeric)
  keep <- logical(length(starts))
  out <- list()
  for (i in seq_along(starts)) {
    s <- starts[[i]]
    duplicate <- any(vapply(out, function(z)
      length(z) == length(s) && max(abs(z - s)) <= tol, logical(1)))
    if (!duplicate) {
      out[[length(out) + 1L]] <- s
      keep[i] <- TRUE
    }
  }
  list(starts = out, keep = keep)
}

.bpss_unique_numeric_tol <- function(x, tol = 1e-12) {
  x <- sort(as.numeric(x[is.finite(x)]))
  if (!length(x)) return(x)
  keep <- c(TRUE, diff(x) > tol * pmax(1, abs(x[-length(x)]), abs(x[-1L])))
  x[keep]
}

.bpss_audit_values <- function(ll_tol, reference_logLik = NA_real_) {
  abs_tol <- as.numeric(ll_tol)[1L]
  rel_tol <- attr(ll_tol, "rel_tol") %||% 1e-8
  rho_tol <- attr(ll_tol, "rho_tol") %||% 0.01
  eta_se_tol <- attr(ll_tol, "eta_se_tol") %||% 0.05
  top_starts <- attr(ll_tol, "top_starts") %||% 3L
  expand <- attr(ll_tol, "expand") %||% TRUE
  effective <- max(abs_tol,
                   rel_tol * max(1, abs(reference_logLik)),
                   na.rm = TRUE)
  list(abs_tol = abs_tol, rel_tol = rel_tol,
       rho_tol = rho_tol, eta_se_tol = eta_se_tol,
       top_starts = as.integer(top_starts), expand = isTRUE(expand),
       effective_ll_tol = effective)
}

# Unrestricted maximum-likelihood fit from a single start.
bpss_fit_unrestricted <- function(D, XS, XO, start = NULL, eta_start = 0,
                                  eps = 1e-12, rcond_tol = 1e-12,
                                  coef_bound = 12,
                                  eta_bound = atanh(0.995),
                                  maxit = 1000, optim_factr = 1e7,
                                  optim_pgtol = 1e-8,
                                  stationarity_tol = 1e-4,
                                  boundary_fraction = 0.99, ...) {
  dd <- theta_dims(XS, XO)
  if (is.null(start))
    start <- make_marginal_start(D, XS, XO, eta_start = eta_start,
                                 coef_bound = min(coef_bound, 8))
  start <- as.numeric(start)
  if (length(start) != dd$ptheta) stop("start has wrong length.")

  ev <- .bpss_make_observed_evaluator(D, XS, XO)
  lower <- c(rep(-coef_bound, dd$ppsi), -eta_bound)
  upper <- c(rep( coef_bound, dd$ppsi),  eta_bound)
  fit <- .fit_bpss_optim(
    start, ev$fn, ev$gr, maxit = maxit,
    factr = optim_factr, pgtol = optim_pgtol,
    stationarity_tol = stationarity_tol,
    lower = lower, upper = upper
  )
  if (!isTRUE(fit$feasible) || length(fit$par) != dd$ptheta)
    return(make_failed_fit(
      "optim_error", fit$error %||% "bounded optimisation failed",
      dd$ptheta, dd$ppsi, theta_attempted = start, optim = fit
    ))

  theta_hat <- fit$par
  names(theta_hat) <- make_theta_names(XS, XO)
  ll <- tryCatch(ev$value(theta_hat)$logLik, error = function(e) NA_real_)
  info <- bpss_info_safe(theta_hat, XS, XO, eps = eps,
                                  rcond_tol = rcond_tol)
  info_ok <- !is.null(info) && !isTRUE(info$probability_floor_hit) &&
    is.finite(info$J) && info$J > 0 &&
    is.finite(info$I_etaeta) && info$I_etaeta > 0
  likelihood_ok <- is.finite(ll)
  optim_ok <- isTRUE(fit$ok)
  bnd <- .bpss_boundary_flags(theta_hat, lower, upper,
                              eta_index = dd$idx_eta,
                              boundary_fraction = boundary_fraction)
  ok <- optim_ok && likelihood_ok && info_ok
  status <- if (!isTRUE(fit$optim_converged)) "optim_nonconvergence"
    else if (!isTRUE(fit$stationary)) "optim_nonstationary"
    else if (!likelihood_ok) "bad_likelihood"
    else if (!info_ok) "bad_information"
    else "ok"

  list(
    ok = ok, optim_ok = optim_ok, likelihood_ok = likelihood_ok,
    info_ok = info_ok, score_ok = NA,
    status = status, message = fit$message %||% NA_character_,
    convergence = fit$convergence, stationary = fit$stationary,
    projected_gradient_norm = fit$projected_gradient_norm,
    projected_gradient_norm_raw = fit$projected_gradient_norm * nrow(D),
    theta = theta_hat, psi = theta_hat[dd$idx_psi],
    logLik = ll, eta_hat = unname(theta_hat[dd$idx_eta]),
    rho_hat = unname(tanh(theta_hat[dd$idx_eta])),
    J = if (info_ok) unname(info$J) else NA_real_,
    I_etaeta = if (info_ok) unname(info$I_etaeta) else NA_real_,
    Delta = if (info_ok) unname(info$Delta) else NA_real_,
    se_eta = if (info_ok) sqrt(1 / info$J) else NA_real_,
    boundary_hit = bnd$boundary_hit,
    eta_boundary_hit = bnd$eta_boundary_hit,
    coef_boundary_hit = bnd$coef_boundary_hit,
    any_at_bound = bnd$any_at_bound,
    eta_at_bound = bnd$eta_at_bound,
    coef_at_bound = bnd$coef_at_bound,
    boundary_components = bnd$hit,
    active_bound_components = bnd$at_bound,
    wald_usable = ok && !bnd$eta_boundary_hit && !bnd$coef_boundary_hit,
    score_usable = NA, info = info, optim = fit, fit_role = "default"
  )
}

# Constrained fit at fixed eta = eta_star; the efficient score is evaluated
# at the fitted nuisance coefficients.
bpss_fit_constrained <- function(eta_star, D, XS, XO,
                                 start_psi = NULL,
                                 eta_start_for_marginal = 0,
                                 eps = 1e-12, rcond_tol = 1e-12,
                                 coef_bound = 12, maxit = 1000,
                                 optim_factr = 1e7,
                                 optim_pgtol = 1e-8,
                                 stationarity_tol = 1e-4,
                                 boundary_fraction = 0.99, ...) {
  dd <- theta_dims(XS, XO)
  if (!is.finite(eta_star)) stop("eta_star must be finite.")
  if (is.null(start_psi)) {
    st <- make_marginal_start(D, XS, XO,
                              eta_start = eta_start_for_marginal,
                              coef_bound = min(coef_bound, 8))
    start_psi <- st[dd$idx_psi]
  }
  start_psi <- as.numeric(start_psi)
  if (length(start_psi) != dd$ppsi) stop("start_psi has wrong length.")

  ev <- .bpss_make_observed_evaluator(D, XS, XO, eta_fixed = eta_star)
  lower <- rep(-coef_bound, dd$ppsi)
  upper <- rep( coef_bound, dd$ppsi)
  fit <- .fit_bpss_optim(
    start_psi, ev$fn, ev$gr, maxit = maxit,
    factr = optim_factr, pgtol = optim_pgtol,
    stationarity_tol = stationarity_tol,
    lower = lower, upper = upper
  )
  if (!isTRUE(fit$feasible) || length(fit$par) != dd$ppsi)
    return(make_failed_fit(
      "constrained_optim_error", fit$error %||% "bounded optimisation failed",
      dd$ptheta, dd$ppsi, psi_attempted = start_psi,
      eta_star = eta_star, optim = fit
    ))

  psi_hat <- fit$par
  theta_hat <- make_theta_from_psi_eta(psi_hat, eta_star, XS, XO)
  ll <- tryCatch(ev$value(psi_hat)$logLik, error = function(e) NA_real_)
  stat_error <- NA_character_
  stat <- tryCatch(
    bpss_efficient_score_stat(theta_hat, XS, XO, D,
                              eps = eps, rcond_tol = rcond_tol),
    error = function(e) {
      stat_error <<- conditionMessage(e)
      NULL
    }
  )
  stat_ok <- !is.null(stat) && is.finite(stat$J) && stat$J > 0 &&
    is.finite(stat$T_score)
  likelihood_ok <- is.finite(ll)
  optim_ok <- isTRUE(fit$ok)
  bnd <- .bpss_boundary_flags(psi_hat, lower, upper,
                              eta_index = NULL,
                              boundary_fraction = boundary_fraction)
  ok <- optim_ok && likelihood_ok && stat_ok
  status <- if (!isTRUE(fit$optim_converged)) "constrained_nonconvergence"
    else if (!isTRUE(fit$stationary)) "constrained_nonstationary"
    else if (!likelihood_ok) "bad_likelihood"
    else if (!stat_ok) "bad_information_or_score"
    else "ok"

  info_obj <- if (!is.null(stat)) stat$info else NULL
  # Smallest fitted cell probability, also when the statistic was not computed.
  min_prob <- if (!is.null(info_obj)) info_obj$min_cell_probability else
    tryCatch(min(bpss_cell_prob(theta_hat, XS, XO)), error = function(e) NA_real_)
  list(
    ok = ok, optim_ok = optim_ok, likelihood_ok = likelihood_ok,
    info_ok = stat_ok, score_ok = stat_ok,
    status = status,
    message = if (!is.na(stat_error)) stat_error else fit$message %||% NA_character_,
    information_error = stat_error,
    convergence = fit$convergence, stationary = fit$stationary,
    projected_gradient_norm = fit$projected_gradient_norm,
    projected_gradient_norm_raw = fit$projected_gradient_norm * nrow(D),
    eta_star = eta_star, theta = theta_hat,
    psi = psi_hat, psi_hat = psi_hat, logLik = ll,
    U = if (!is.null(stat)) stat$U else NA_real_,
    J = if (!is.null(stat)) stat$J else NA_real_,
    J_schur = if (!is.null(info_obj)) info_obj$J_schur else NA_real_,
    I_etaeta = if (!is.null(stat)) stat$I_etaeta else NA_real_,
    T_score = if (!is.null(stat)) stat$T_score else NA_real_,
    Delta = if (!is.null(stat)) stat$Delta else NA_real_,
    rcond_psipsi = if (!is.null(info_obj)) info_obj$rcond_psipsi else NA_real_,
    min_cell_probability = min_prob,
    probability_floor_hit = if (!is.null(info_obj))
      info_obj$probability_floor_hit else if (is.finite(min_prob))
      min_prob < eps else NA,
    boundary_hit = bnd$boundary_hit,
    eta_boundary_hit = FALSE,
    coef_boundary_hit = bnd$coef_boundary_hit,
    any_at_bound = bnd$any_at_bound,
    eta_at_bound = FALSE,
    coef_at_bound = bnd$coef_at_bound,
    boundary_components = bnd$hit,
    active_bound_components = bnd$at_bound,
    score_usable = ok && !bnd$coef_boundary_hit,
    wald_usable = NA, stat = stat, optim = fit,
    fit_role = "constrained"
  )
}

# Multi-start audit at fixed eta. Parameter differences describe possible
# nuisance nonuniqueness, but a continuation fit is called discordant only when
# an independent start produces a materially larger likelihood.
bpss_fit_constrained_audited <- function(
    eta_star, D, XS, XO, start_psi = NULL, extra_starts = NULL,
    eps = 1e-12, rcond_tol = 1e-12, coef_bound = 12,
    maxit = 1000, optim_factr = 1e7, optim_pgtol = 1e-8,
    stationarity_tol = 1e-4, boundary_fraction = 0.99,
    ll_tol = structure(1e-3, rel_tol = 1e-8),
    par_tol = 1e-3, ...) {
  dd <- theta_dims(XS, XO)
  marginal <- make_marginal_start(D, XS, XO, eta_start = eta_star,
                                  coef_bound = min(coef_bound, 8))[dd$idx_psi]
  starts <- list(); labels <- character(0)
  add_start <- function(x, label) {
    if (!is.null(x) && length(x) == dd$ppsi && all(is.finite(x))) {
      starts[[length(starts) + 1L]] <<- as.numeric(x)
      labels <<- c(labels, label)
    }
  }
  add_start(start_psi, "continuation")
  add_start(marginal, "marginal")
  if (!is.null(extra_starts)) {
    if (is.numeric(extra_starts)) extra_starts <- list(extra_starts)
    for (j in seq_along(extra_starts))
      add_start(extra_starts[[j]], paste0("extra_", j))
  }
  uu <- .bpss_unique_starts(starts)
  labels <- labels[uu$keep]; starts <- uu$starts
  fits <- lapply(starts, function(st)
    bpss_fit_constrained(
      eta_star, D, XS, XO, start_psi = st, eps = eps,
      rcond_tol = rcond_tol, coef_bound = coef_bound, maxit = maxit,
      optim_factr = optim_factr, optim_pgtol = optim_pgtol,
      stationarity_tol = stationarity_tol,
      boundary_fraction = boundary_fraction
    ))
  eligible <- vapply(fits, function(z)
    isTRUE(z$optim_ok) && isTRUE(z$likelihood_ok) && is.finite(z$logLik),
    logical(1))
  if (!any(eligible)) {
    ans <- fits[[1L]]
    ans$audit_ok <- FALSE
    ans$audit <- list(labels = labels, fits = fits, eligible = eligible,
                      best = NA_integer_, branch_ambiguous = NA,
                      start_discordance = NA)
    return(ans)
  }
  idx <- which(eligible)
  ll <- vapply(fits[idx], `[[`, numeric(1), "logLik")
  best <- idx[which.max(ll)]
  ans <- fits[[best]]
  av <- .bpss_audit_values(ll_tol, ans$logLik)
  near <- idx[(max(ll) - ll) <= av$effective_ll_tol]
  branch_ambiguous <- FALSE
  if (length(near) > 1L) {
    ref <- fits[[near[1L]]]$psi
    branch_ambiguous <- any(vapply(near[-1L], function(k)
      max(abs(fits[[k]]$psi - ref)) > par_tol, logical(1)))
  }
  continuation_idx <- match("continuation", labels)
  ll_gain <- if (is.finite(continuation_idx) && eligible[continuation_idx])
    ans$logLik - fits[[continuation_idx]]$logLik else NA_real_
  start_discordance <- is.finite(ll_gain) &&
    ll_gain > av$effective_ll_tol

  ans$audit_ok <- isTRUE(ans$optim_ok) && isTRUE(ans$likelihood_ok)
  ans$branch_ambiguous <- branch_ambiguous
  ans$start_discordance <- start_discordance
  ans$audit <- list(
    labels = labels, fits = fits, eligible = eligible, best = best,
    best_label = labels[best], best_logLik = ans$logLik,
    logLik_range = range(ll), near_best = near,
    branch_ambiguous = branch_ambiguous,
    start_discordance = start_discordance,
    continuation_logLik_gain = ll_gain,
    effective_ll_tol = av$effective_ll_tol,
    ll_tol = av$abs_tol, ll_rel_tol = av$rel_tol, par_tol = par_tol
  )
  ans$fit_role <- "constrained_audited"
  ans
}

# Adaptive unrestricted audit. Seven fixed-eta profile fits are retained, but
# unrestricted optimisation is launched first only from the default, the best
# profile starts and the two profile extremes. Remaining starts are expanded
# automatically if the initial audit detects a materially different solution or
# a constrained profile value above the unrestricted reference.
bpss_fit_unrestricted_audited <- function(
    D, XS, XO, default_fit = NULL, start = NULL,
    eta_starts = atanh(c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)),
    extra_starts = NULL, use_constrained_starts = TRUE,
    eps = 1e-12, rcond_tol = 1e-12, coef_bound = 12,
    eta_bound = atanh(0.995), maxit = 1000,
    optim_factr = 1e7, optim_pgtol = 1e-8,
    stationarity_tol = 1e-4, boundary_fraction = 0.99,
    ll_tol = structure(1e-3, rel_tol = 1e-8, rho_tol = 0.01,
                       eta_se_tol = 0.05, top_starts = 3L, expand = TRUE),
    par_tol = 1e-3, ...) {
  dd <- theta_dims(XS, XO)
  if (is.null(default_fit)) {
    default_fit <- bpss_fit_unrestricted(
      D, XS, XO, start = start, eps = eps, rcond_tol = rcond_tol,
      coef_bound = coef_bound, eta_bound = eta_bound, maxit = maxit,
      optim_factr = optim_factr, optim_pgtol = optim_pgtol,
      stationarity_tol = stationarity_tol,
      boundary_fraction = boundary_fraction
    )
  }
  marginal <- make_marginal_start(D, XS, XO, eta_start = 0,
                                  coef_bound = min(coef_bound, 8))
  eta_starts <- sort(.bpss_unique_numeric_tol(
    pmin(pmax(as.numeric(eta_starts), -eta_bound), eta_bound),
    tol = 1e-12
  ))

  constrained_checks <- vector("list", length(eta_starts))
  profile_starts <- vector("list", length(eta_starts))
  for (j in seq_along(eta_starts)) {
    e <- eta_starts[j]
    if (isTRUE(use_constrained_starts)) {
      fc <- bpss_fit_constrained(
        e, D, XS, XO, start_psi = marginal[dd$idx_psi],
        eps = eps, rcond_tol = rcond_tol, coef_bound = coef_bound,
        maxit = maxit, optim_factr = optim_factr,
        optim_pgtol = optim_pgtol,
        stationarity_tol = stationarity_tol,
        boundary_fraction = boundary_fraction
      )
      constrained_checks[[j]] <- fc
      if (isTRUE(fc$likelihood_ok) && all(is.finite(fc$theta)))
        profile_starts[[j]] <- fc$theta
      else
        profile_starts[[j]] <- make_theta_from_psi_eta(
          marginal[dd$idx_psi], e, XS, XO)
    } else {
      profile_starts[[j]] <- make_theta_from_psi_eta(
        marginal[dd$idx_psi], e, XS, XO)
    }
  }

  profile_ll <- vapply(constrained_checks, function(z)
    if (!is.null(z) && isTRUE(z$likelihood_ok)) z$logLik else -Inf,
    numeric(1))
  av0 <- .bpss_audit_values(ll_tol,
                            if (is.finite(default_fit$logLik))
                              default_fit$logLik else 0)
  top <- head(order(profile_ll, decreasing = TRUE), av0$top_starts)
  initial_profile <- unique(c(top, 1L, length(profile_starts)))

  starts <- list(); labels <- character(0)
  add_start <- function(x, label) {
    if (!is.null(x) && length(x) == dd$ptheta && all(is.finite(x))) {
      starts[[length(starts) + 1L]] <<- as.numeric(x)
      labels <<- c(labels, label)
    }
  }
  if (!is.null(default_fit$theta) && all(is.finite(default_fit$theta)))
    add_start(default_fit$theta, "default")
  add_start(start, "user")
  for (j in initial_profile)
    add_start(profile_starts[[j]], paste0("profile_", j))
  if (!is.null(extra_starts)) {
    if (is.numeric(extra_starts)) extra_starts <- list(extra_starts)
    for (j in seq_along(extra_starts))
      add_start(extra_starts[[j]], paste0("extra_", j))
  }
  uu <- .bpss_unique_starts(starts)
  labels <- labels[uu$keep]; starts <- uu$starts

  fit_from_start <- function(st, label) {
    if (identical(label, "default") &&
        !is.null(default_fit$theta) &&
        max(abs(st - default_fit$theta)) <= 1e-10) return(default_fit)
    bpss_fit_unrestricted(
      D, XS, XO, start = st, eps = eps, rcond_tol = rcond_tol,
      coef_bound = coef_bound, eta_bound = eta_bound, maxit = maxit,
      optim_factr = optim_factr, optim_pgtol = optim_pgtol,
      stationarity_tol = stationarity_tol,
      boundary_fraction = boundary_fraction
    )
  }
  fits <- lapply(seq_along(starts), function(j)
    fit_from_start(starts[[j]], labels[j]))

  eligible_now <- function(ff) vapply(ff, function(z)
    isTRUE(z$optim_ok) && isTRUE(z$likelihood_ok) && is.finite(z$logLik),
    logical(1))
  eligible <- eligible_now(fits)
  if (!any(eligible)) {
    ans <- default_fit
    ans$audit_ok <- FALSE
    ans$reference_dominated <- NA
    ans$default_discordance <- NA
    ans$default_logLik_gain <- NA_real_
    ans$default_parameter_gap <- NA_real_
    ans$multiple_near_maxima <- NA
    ans$default_fit <- default_fit
    ans$audit <- list(labels = labels, starts = starts, fits = fits,
                      eligible = eligible, best = NA_integer_)
    ans$fit_role <- "audited_failed"
    return(ans)
  }

  initial_ll <- vapply(fits[eligible], `[[`, numeric(1), "logLik")
  initial_rho <- vapply(fits[eligible], `[[`, numeric(1), "rho_hat")
  initial_best_ll <- max(initial_ll)
  avi <- .bpss_audit_values(ll_tol, initial_best_ll)
  max_constrained_ll <- if (any(is.finite(profile_ll)))
    max(profile_ll[is.finite(profile_ll)]) else NA_real_
  need_expand <- isTRUE(avi$expand) && (
    (is.finite(max_constrained_ll) &&
       max_constrained_ll > initial_best_ll + avi$effective_ll_tol) ||
    diff(range(initial_rho)) > avi$rho_tol ||
    diff(range(initial_ll)) > avi$effective_ll_tol
  )

  if (need_expand) {
    remaining <- setdiff(seq_along(profile_starts), initial_profile)
    for (j in remaining) {
      st <- profile_starts[[j]]
      duplicate <- any(vapply(starts, function(z)
        max(abs(z - st)) <= 1e-10, logical(1)))
      if (!duplicate) {
        starts[[length(starts) + 1L]] <- st
        labels <- c(labels, paste0("profile_", j))
        fits[[length(fits) + 1L]] <- fit_from_start(st, tail(labels, 1L))
      }
    }
  }

  eligible <- eligible_now(fits)
  idx <- which(eligible)
  ll <- vapply(fits[idx], `[[`, numeric(1), "logLik")
  best <- idx[which.max(ll)]
  ans <- fits[[best]]
  av <- .bpss_audit_values(ll_tol, ans$logLik)
  near <- idx[(max(ll) - ll) <= av$effective_ll_tol]
  near_rho <- vapply(fits[near], `[[`, numeric(1), "rho_hat")
  near_eta <- vapply(fits[near], `[[`, numeric(1), "eta_hat")
  se_eta <- ans$se_eta
  multiple_near_maxima <- length(near) > 1L && (
    diff(range(near_rho)) > av$rho_tol ||
    (is.finite(se_eta) && se_eta > 0 &&
       diff(range(near_eta)) / se_eta > av$eta_se_tol)
  )

  max_constrained_ll <- if (any(is.finite(profile_ll)))
    max(profile_ll[is.finite(profile_ll)]) else NA_real_
  reference_dominated <- is.finite(max_constrained_ll) &&
    max_constrained_ll > ans$logLik + av$effective_ll_tol

  default_ll <- if (isTRUE(default_fit$likelihood_ok))
    default_fit$logLik else NA_real_
  ll_gain <- if (is.finite(default_ll)) ans$logLik - default_ll else NA_real_
  parameter_gap <- if (!is.null(default_fit$theta) &&
                       all(is.finite(default_fit$theta)))
    max(abs(ans$theta - default_fit$theta)) else NA_real_
  rho_gap <- if (is.finite(default_fit$rho_hat))
    abs(ans$rho_hat - default_fit$rho_hat) else NA_real_
  eta_gap <- if (is.finite(default_fit$eta_hat))
    abs(ans$eta_hat - default_fit$eta_hat) else NA_real_
  eta_gap_se <- if (is.finite(eta_gap) && is.finite(se_eta) && se_eta > 0)
    eta_gap / se_eta else NA_real_
  ll_material <- is.finite(ll_gain) && ll_gain > av$effective_ll_tol
  rho_material <- is.finite(rho_gap) && rho_gap > av$rho_tol
  eta_material <- is.finite(eta_gap_se) && eta_gap_se > av$eta_se_tol
  default_discordance <- ll_material || rho_material || eta_material
  numerical_difference <- (is.finite(ll_gain) && ll_gain > 1e-7) ||
    (is.finite(parameter_gap) && parameter_gap > 1e-4)

  ans$audit_ok <- isTRUE(ans$optim_ok) && isTRUE(ans$likelihood_ok) &&
    !reference_dominated
  ans$reference_dominated <- reference_dominated
  ans$default_discordance <- default_discordance
  ans$default_numerical_difference <- numerical_difference
  ans$default_logLik_gain <- ll_gain
  ans$default_parameter_gap <- parameter_gap
  ans$default_rho_gap <- rho_gap
  ans$default_eta_gap <- eta_gap
  ans$default_eta_gap_se <- eta_gap_se
  ans$multiple_near_maxima <- multiple_near_maxima
  ans$default_fit <- default_fit
  ans$audit <- list(
    labels = labels, starts = starts, fits = fits,
    eligible = eligible, best = best, best_label = labels[best],
    best_logLik = ans$logLik, logLik_range = range(ll),
    near_best = near, multiple_near_maxima = multiple_near_maxima,
    constrained_checks = constrained_checks,
    max_constrained_logLik = max_constrained_ll,
    reference_dominated = reference_dominated,
    default_logLik = default_ll, default_logLik_gain = ll_gain,
    default_parameter_gap = parameter_gap,
    default_rho_gap = rho_gap, default_eta_gap = eta_gap,
    default_eta_gap_se = eta_gap_se,
    default_discordance = default_discordance,
    default_numerical_difference = numerical_difference,
    effective_ll_tol = av$effective_ll_tol,
    ll_tol = av$abs_tol, ll_rel_tol = av$rel_tol,
    rho_tol = av$rho_tol, eta_se_tol = av$eta_se_tol,
    par_tol = par_tol, expanded = need_expand
  )
  ans$fit_role <- "audited"
  ans
}

# Probit fit by the same bounded optimiser.
.bpss_fit_probit_bounded <- function(X, y, start = NULL, coef_bound = 12,
                                      maxit = 1000, optim_factr = 1e7,
                                      optim_pgtol = 1e-8,
                                      stationarity_tol = 1e-4,
                                      boundary_fraction = 0.99) {
  X <- as.matrix(X); y <- as.numeric(y)
  if (nrow(X) != length(y)) stop("X and y have incompatible dimensions.")
  if (is.null(start))
    start <- safe_glm_probit_coef(X, y, coef_bound = min(coef_bound, 8))
  start <- as.numeric(start)
  if (length(start) != ncol(X)) stop("probit start has wrong length.")

  n_obs <- nrow(X)
  objective <- function(b) {
    z <- drop(X %*% b)
    p <- stats::pnorm(z)
    q <- ifelse(y == 1, p, 1 - p)
    if (any(!is.finite(q)) || any(q <= 0))
      return(.Machine$double.xmax / 100)
    -sum(log(q)) / n_obs
  }
  gradient <- function(b) {
    z <- drop(X %*% b)
    p <- stats::pnorm(z)
    phi <- stats::dnorm(z)
    q <- ifelse(y == 1, p, 1 - p)
    if (any(!is.finite(q)) || any(q <= 0))
      stop("zero probability in bounded probit gradient.")
    s <- ifelse(y == 1, phi / p, -phi / (1 - p))
    -drop(crossprod(X, s)) / n_obs
  }
  lower <- rep(-coef_bound, ncol(X))
  upper <- rep( coef_bound, ncol(X))
  fit <- .fit_bpss_optim(
    start, objective, gradient, maxit = maxit,
    factr = optim_factr, pgtol = optim_pgtol,
    stationarity_tol = stationarity_tol,
    lower = lower, upper = upper
  )
  bnd <- .bpss_boundary_flags(fit$par, lower, upper,
                              boundary_fraction = boundary_fraction)
  list(ok = isTRUE(fit$ok), par = fit$par,
       logLik = if (is.finite(fit$value))
         -fit$value * n_obs else NA_real_,
       optim = fit, boundary_hit = bnd$boundary_hit,
       any_at_bound = bnd$any_at_bound,
       coef_at_bound = bnd$coef_at_bound,
       projected_gradient_norm = fit$projected_gradient_norm,
       projected_gradient_norm_raw =
         fit$projected_gradient_norm * n_obs,
       stationary = fit$stationary)
}

# Constrained fit at eta = 0 by factorisation.
bpss_fit_null_factorized <- function(
    D, XS, XO, eps = 1e-12, rcond_tol = 1e-12,
    coef_bound = 12, maxit = 1000, optim_factr = 1e7,
    optim_pgtol = 1e-8, stationarity_tol = 1e-4,
    boundary_fraction = 0.99, check_generic = TRUE,
    generic_ll_tol = 1e-3, generic_par_tol = 1e-3, ...) {
  D <- as.matrix(D); XS <- as.matrix(XS); XO <- as.matrix(XO)
  dd <- theta_dims(XS, XO)
  dy <- D_to_YS_YO(D)
  S <- dy$YS; Y <- dy$YO
  selected <- S == 1
  if (sum(selected) < ncol(XO) + 1L)
    return(make_failed_fit("null_too_few_selected",
                           "too few selected observations for the outcome probit",
                           dd$ptheta, dd$ppsi, eta_star = 0))

  gamma_fit <- .bpss_fit_probit_bounded(
    XS, S, coef_bound = coef_bound, maxit = maxit,
    optim_factr = optim_factr, optim_pgtol = optim_pgtol,
    stationarity_tol = stationarity_tol,
    boundary_fraction = boundary_fraction
  )
  beta_fit <- .bpss_fit_probit_bounded(
    XO[selected, , drop = FALSE], Y[selected],
    coef_bound = coef_bound, maxit = maxit,
    optim_factr = optim_factr, optim_pgtol = optim_pgtol,
    stationarity_tol = stationarity_tol,
    boundary_fraction = boundary_fraction
  )
  theta <- c(gamma_fit$par, beta_fit$par, eta = 0)
  names(theta) <- make_theta_names(XS, XO)
  ll_components <- gamma_fit$logLik + beta_fit$logLik
  ll_full <- bpss_loglik(theta, XS, XO, D = D, eps = eps)
  optim_ok <- isTRUE(gamma_fit$ok) && isTRUE(beta_fit$ok)
  likelihood_ok <- is.finite(ll_full) && is.finite(ll_components)
  stat <- tryCatch(
    bpss_efficient_score_stat(theta, XS, XO, D,
                              eps = eps, rcond_tol = rcond_tol),
    error = function(e) NULL
  )
  stat_ok <- !is.null(stat) && is.finite(stat$J) && stat$J > 0 &&
             is.finite(stat$T_score)
  coef_boundary_hit <- isTRUE(gamma_fit$boundary_hit) ||
                       isTRUE(beta_fit$boundary_hit)
  coef_at_bound <- isTRUE(gamma_fit$coef_at_bound) ||
                   isTRUE(beta_fit$coef_at_bound)

  generic <- NULL
  generic_ll_difference <- generic_parameter_difference <- NA_real_
  generic_agrees <- NA
  if (isTRUE(check_generic)) {
    generic <- bpss_fit_constrained(
      0, D, XS, XO, start_psi = theta[dd$idx_psi],
      eps = eps, rcond_tol = rcond_tol, coef_bound = coef_bound,
      maxit = maxit, optim_factr = optim_factr,
      optim_pgtol = optim_pgtol,
      stationarity_tol = stationarity_tol,
      boundary_fraction = boundary_fraction
    )
    if (isTRUE(generic$optim_ok) && isTRUE(generic$likelihood_ok)) {
      generic_ll_difference <- generic$logLik - ll_full
      generic_parameter_difference <- max(abs(generic$psi - theta[dd$idx_psi]))
      generic_agrees <- abs(generic_ll_difference) <= generic_ll_tol &&
                        generic_parameter_difference <= generic_par_tol
    } else {
      generic_agrees <- FALSE
    }
  }

  # At eta = 0 the constrained likelihood factorises, so the two probit fits
  # give the constrained maximum. The generic constrained fit is a numerical
  # cross-check only: its outcome is reported in generic_agrees and does not
  # change ok or status.
  ok <- optim_ok && likelihood_ok && stat_ok
  status <- if (!optim_ok) "null_nonconvergence"
            else if (!likelihood_ok) "null_bad_likelihood"
            else if (!stat_ok) "null_bad_information_or_score"
            else "ok"

  list(
    ok = ok, optim_ok = optim_ok, likelihood_ok = likelihood_ok,
    info_ok = stat_ok, score_ok = stat_ok,
    status = status, message = NA_character_,
    theta = theta, psi = theta[dd$idx_psi], psi_hat = theta[dd$idx_psi],
    eta_star = 0, eta_hat = 0, rho_hat = 0,
    logLik = ll_full, component_logLik = ll_components,
    factorization_error = ll_full - ll_components,
    U = if (!is.null(stat)) stat$U else NA_real_,
    J = if (!is.null(stat)) stat$J else NA_real_,
    I_etaeta = if (!is.null(stat)) stat$I_etaeta else NA_real_,
    T_score = if (!is.null(stat)) stat$T_score else NA_real_,
    Delta = if (!is.null(stat)) stat$Delta else NA_real_,
    boundary_hit = coef_boundary_hit,
    eta_boundary_hit = FALSE,
    coef_boundary_hit = coef_boundary_hit,
    any_at_bound = coef_at_bound,
    eta_at_bound = FALSE,
    coef_at_bound = coef_at_bound,
    score_usable = ok && !coef_boundary_hit,
    wald_usable = NA,
    stationary = optim_ok,
    projected_gradient_norm = max(gamma_fit$projected_gradient_norm,
                                  beta_fit$projected_gradient_norm),
    gamma_fit = gamma_fit, beta_fit = beta_fit,
    generic = generic, generic_agrees = generic_agrees,
    generic_logLik_difference = generic_ll_difference,
    generic_parameter_difference = generic_parameter_difference,
    stat = stat, fit_role = "null_factorized"
  )
}
