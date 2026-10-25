# Inference on eta = atanh(rho): grid inversion of the constrained score
# statistic with refinement at changes in the decision, the likelihood-ratio
# profile on the same grid, the test of ignorable selection (rho = 0), the
# score test at a specified correlation, and the Wald interval.

interval_total_length <- function(intervals) {
  if (is.null(intervals) || nrow(intervals) == 0L) return(0)
  sum(pmax(0, intervals$upper - intervals$lower))
}

transform_eta_intervals_to_rho <- function(intervals) {
  if (is.null(intervals) || nrow(intervals) == 0L)
    return(data.frame(lower = numeric(0), upper = numeric(0)))
  data.frame(lower = tanh(intervals$lower), upper = tanh(intervals$upper))
}

grid_intervals_with_edges <- function(grid, accept) {
  grid <- as.numeric(grid); accept <- as.logical(accept)
  ok <- which(accept & is.finite(grid))
  if (length(ok) == 0L)
    return(list(intervals = data.frame(lower = numeric(0), upper = numeric(0)),
                index_runs = list(), touches_left = FALSE, touches_right = FALSE, touches_edge = FALSE))
  runs <- split(ok, cumsum(c(TRUE, diff(ok) != 1L)))
  intervals <- do.call(rbind, lapply(runs, function(ii) data.frame(lower = min(grid[ii]), upper = max(grid[ii]))))
  rownames(intervals) <- NULL
  touches_left  <- any(vapply(runs, function(ii) min(ii) == 1L, logical(1)))
  touches_right <- any(vapply(runs, function(ii) max(ii) == length(grid), logical(1)))
  list(intervals = intervals, index_runs = runs, touches_left = touches_left,
       touches_right = touches_right, touches_edge = touches_left || touches_right)
}

wald_eta_interval <- function(eta_hat, J, level = 0.95) {
  if (!is.finite(eta_hat) || !is.finite(J) || J <= 0)
    return(data.frame(lower = numeric(0), upper = numeric(0)))
  z <- stats::qnorm(1 - (1 - level) / 2); half <- z / sqrt(J)
  data.frame(lower = eta_hat - half, upper = eta_hat + half)
}

# Refine one accepted run [lo, hi] of grid indices to continuous endpoints by
# root-finding on margin_fun(eta, warm_psi) (nonnegative = not rejected), when
# the neighbouring grid value is resolved and rejected.
.bpss_refine_run <- function(lo, hi, eta_grid, accept, resolved, psi_list,
                             margin_fun, tol = 1e-6) {
  ngrid <- length(eta_grid)
  lower <- eta_grid[lo]; upper <- eta_grid[hi]
  edge_left <- (lo == 1L); edge_right <- (hi == ngrid)
  # left endpoint: bracket [eta[lo-1], eta[lo]] if lo>1 and lo-1 resolved & rejecting
  if (lo > 1L && isTRUE(resolved[lo - 1L]) && !isTRUE(accept[lo - 1L])) {
    warm <- psi_list[[lo]]
    root <- tryCatch(stats::uniroot(function(e) margin_fun(e, warm),
                                    lower = eta_grid[lo - 1L], upper = eta_grid[lo],
                                    tol = tol)$root, error = function(e) NA_real_)
    if (is.finite(root)) lower <- root
  }
  # right endpoint: bracket [eta[hi], eta[hi+1]] if hi<ngrid and hi+1 resolved & rejecting
  if (hi < ngrid && isTRUE(resolved[hi + 1L]) && !isTRUE(accept[hi + 1L])) {
    warm <- psi_list[[hi]]
    root <- tryCatch(stats::uniroot(function(e) margin_fun(e, warm),
                                    lower = eta_grid[hi], upper = eta_grid[hi + 1L],
                                    tol = tol)$root, error = function(e) NA_real_)
    if (is.finite(root)) upper <- root
  }
  c(lower = lower, upper = upper, edge_left = edge_left, edge_right = edge_right)
}

.bpss_grid_audit_spec <- function(x) {
  if (isTRUE(x)) return(list(mode = "all", stride = 1L))
  if (identical(x, FALSE) || is.null(x))
    return(list(mode = "none", stride = Inf))
  if (is.list(x)) {
    mode <- match.arg(x$mode %||% "adaptive", c("adaptive", "all", "none"))
    stride <- as.integer(x$stride %||% 20L)
    return(list(mode = mode, stride = max(1L, stride)))
  }
  mode <- match.arg(as.character(x)[1L], c("adaptive", "all", "none"))
  stride <- as.integer(attr(x, "stride") %||% 20L)
  list(mode = mode, stride = max(1L, stride))
}

.bpss_unresolved_ranges <- function(tab) {
  ii <- which(!tab$resolved)
  if (!length(ii))
    return(data.frame(lower = numeric(0), upper = numeric(0)))
  runs <- split(ii, cumsum(c(TRUE, diff(ii) != 1L)))
  do.call(rbind, lapply(runs, function(j)
    data.frame(lower = min(tab$rho[j]), upper = max(tab$rho[j]))))
}

# Grid inversion. Constrained fits are obtained by bidirectional continuation
# from the grid value nearest the unrestricted estimate. With the default
# adaptive audit, fits from the separate-probit and unrestricted nuisance starts
# are compared at the two grid ends, the continuation centre, eta = 0, the
# values in audit_eta_values, every stride-th grid value, changes in the
# decision, values at which the statistic is not computed and their
# neighbours, and values whose constrained log-likelihood exceeds the
# unrestricted reference. constrained_audit = "all" audits every value and
# "none" uses continuation only.
bpss_grid_invert <- function(
    D, XS, XO, eta_grid, unrestricted = NULL,
    profile_reference = NULL, level = 0.95,
    eps = 1e-12, rcond_tol = 1e-12, coef_bound = 12,
    maxit = 1000, optim_factr = 1e7, optim_pgtol = 1e-8,
    stationarity_tol = 1e-4, boundary_fraction = 0.99,
    refine = TRUE, refine_tol = 1e-6,
    keep_fits = TRUE, constrained_audit = "adaptive",
    audit_ll_tol = structure(1e-3, rel_tol = 1e-8),
    audit_par_tol = 1e-3, audit_eta_values = numeric(0),
    require_audited_profile = TRUE, profile_lr_tol = 1e-3) {
  eta_grid <- .bpss_unique_numeric_tol(eta_grid, tol = 1e-12)
  ngrid <- length(eta_grid)
  dd <- theta_dims(XS, XO)
  alpha <- 1 - level
  cutoff <- stats::qchisq(level, df = 1)
  audit_spec <- .bpss_grid_audit_spec(constrained_audit)
  audit_eta_values <- as.numeric(audit_eta_values)
  audit_eta_values <- audit_eta_values[is.finite(audit_eta_values)]

  marg <- make_marginal_start(D, XS, XO)[dd$idx_psi]
  if (!is.null(unrestricted) && isTRUE(unrestricted$optim_ok) &&
      isTRUE(unrestricted$likelihood_ok)) {
    center_psi <- unrestricted$psi
    center_eta <- unrestricted$eta_hat
  } else {
    center_psi <- marg
    center_eta <- stats::median(eta_grid)
  }

  if (is.null(profile_reference) && !is.null(unrestricted) &&
      isTRUE(unrestricted$audit_ok)) profile_reference <- unrestricted
  profile_requested <- !is.null(profile_reference)
  profile_reference_ok <- profile_requested &&
    isTRUE(profile_reference$optim_ok) &&
    isTRUE(profile_reference$likelihood_ok) &&
    is.finite(profile_reference$logLik) &&
    (!isTRUE(require_audited_profile) || isTRUE(profile_reference$audit_ok))
  ll_ref <- if (profile_reference_ok) profile_reference$logLik else NA_real_

  tab <- data.frame(
    eta = eta_grid, rho = tanh(eta_grid), status = NA_character_,
    ok = FALSE, resolved = FALSE, score_usable = FALSE,
    logLik = NA_real_, U = NA_real_, J = NA_real_, J_schur = NA_real_,
    I_etaeta = NA_real_, Delta = NA_real_, T = NA_real_,
    p_chisq = NA_real_,
    rcond_psipsi = NA_real_, min_cell_probability = NA_real_,
    probability_floor_hit = NA, information_error = NA_character_,
    coef_boundary_hit = NA, projected_gradient_norm = NA_real_,
    audit_performed = FALSE, branch_ambiguous = NA,
    start_discordance = NA, LR_raw = NA_real_, LR = NA_real_,
    profile_resolved = FALSE, reference_dominated = FALSE,
    stringsAsFactors = FALSE
  )
  fits <- vector("list", ngrid)
  psi_list <- vector("list", ngrid)
  fit_cache <- new.env(hash = TRUE, parent = emptyenv())
  cache_key <- function(eta) sprintf("%.14g", eta)
  cache_put <- function(eta, fit) {
    key <- cache_key(eta)
    if (!exists(key, fit_cache, inherits = FALSE)) {
      assign(key, fit, fit_cache)
    } else {
      old <- get(key, fit_cache, inherits = FALSE)
      if (is.finite(fit$logLik) &&
          (!is.finite(old$logLik) || fit$logLik > old$logLik))
        assign(key, fit, fit_cache)
    }
    invisible(fit)
  }
  cache_get <- function(eta) {
    key <- cache_key(eta)
    if (exists(key, fit_cache, inherits = FALSE))
      get(key, fit_cache, inherits = FALSE) else NULL
  }

  fit_single <- function(eta, warm, use_cache = TRUE) {
    if (isTRUE(use_cache)) {
      got <- cache_get(eta)
      if (!is.null(got)) return(got)
    }
    ans <- tryCatch(
      bpss_fit_constrained(
        eta, D, XS, XO, start_psi = warm,
        eps = eps, rcond_tol = rcond_tol,
        coef_bound = coef_bound, maxit = maxit,
        optim_factr = optim_factr, optim_pgtol = optim_pgtol,
        stationarity_tol = stationarity_tol,
        boundary_fraction = boundary_fraction
      ),
      error = function(e) make_failed_fit(
        "constrained_error", conditionMessage(e),
        dd$ptheta, dd$ppsi, eta_star = eta)
    )
    cache_put(eta, ans)
    ans
  }

  update_row <- function(j, fc, audited = FALSE) {
    fits[[j]] <<- fc
    tab$status[j] <<- fc$status
    tab$ok[j] <<- isTRUE(fc$ok)
    resolved <- isTRUE(fc$score_usable) && is.finite(fc$J) && fc$J > 0 &&
      is.finite(fc$T_score)
    tab$resolved[j] <<- resolved
    tab$score_usable[j] <<- isTRUE(fc$score_usable)
    tab$logLik[j] <<- fc$logLik
    tab$U[j] <<- fc$U
    tab$J[j] <<- fc$J
    tab$J_schur[j] <<- fc$J_schur %||% NA_real_
    tab$I_etaeta[j] <<- fc$I_etaeta
    tab$Delta[j] <<- fc$Delta
    tab$T[j] <<- fc$T_score
    tab$rcond_psipsi[j] <<- fc$rcond_psipsi %||% NA_real_
    tab$min_cell_probability[j] <<- fc$min_cell_probability %||% NA_real_
    tab$probability_floor_hit[j] <<- fc$probability_floor_hit %||% NA
    tab$information_error[j] <<- fc$information_error %||% NA_character_
    tab$coef_boundary_hit[j] <<- fc$coef_boundary_hit
    tab$projected_gradient_norm[j] <<- fc$projected_gradient_norm
    tab$audit_performed[j] <<- audited
    tab$branch_ambiguous[j] <<- isTRUE(fc$branch_ambiguous)
    tab$start_discordance[j] <<- isTRUE(fc$start_discordance)
    if (isTRUE(fc$optim_ok) && isTRUE(fc$likelihood_ok))
      psi_list[[j]] <<- fc$psi_hat

    tab$p_chisq[j] <<- NA_real_
    if (resolved)
      tab$p_chisq[j] <<- stats::pchisq(fc$T_score, df = 1,
                                       lower.tail = FALSE)

    # LR = 2 {l(reference) - l_p(eta)} = -2 log R_n(eta).
    tab$LR_raw[j] <<- tab$LR[j] <<- NA_real_
    tab$profile_resolved[j] <<- FALSE
    tab$reference_dominated[j] <<- FALSE
    if (profile_reference_ok && isTRUE(fc$optim_ok) &&
        isTRUE(fc$likelihood_ok)) {
      raw <- 2 * (ll_ref - fc$logLik)
      tab$LR_raw[j] <<- raw
      dominated <- is.finite(raw) && raw < -profile_lr_tol
      tab$reference_dominated[j] <<- dominated
      if (is.finite(raw) && !dominated) {
        tab$LR[j] <<- max(0, raw)
        tab$profile_resolved[j] <<- TRUE
      }
    }
    invisible(fc)
  }

  audit_one <- function(j) {
    base <- fits[[j]]
    extras <- list(marg, center_psi)
    extra_fits <- list()
    for (st in extras) {
      if (!is.null(base$psi_hat) && all(is.finite(base$psi_hat)) &&
          max(abs(st - base$psi_hat)) <= 1e-10) next
      extra_fits[[length(extra_fits) + 1L]] <- fit_single(
        eta_grid[j], st, use_cache = FALSE)
    }
    all_fits <- c(list(base), extra_fits)
    eligible <- vapply(all_fits, function(z)
      isTRUE(z$optim_ok) && isTRUE(z$likelihood_ok) && is.finite(z$logLik),
      logical(1))
    if (!any(eligible)) {
      update_row(j, base, audited = TRUE)
      return(invisible(base))
    }
    idx <- which(eligible)
    ll <- vapply(all_fits[idx], `[[`, numeric(1), "logLik")
    best <- idx[which.max(ll)]
    ans <- all_fits[[best]]
    avj <- .bpss_audit_values(audit_ll_tol, ans$logLik)
    near <- idx[(max(ll) - ll) <= avj$effective_ll_tol]
    branch <- FALSE
    if (length(near) > 1L) {
      ref <- all_fits[[near[1L]]]$psi
      branch <- any(vapply(near[-1L], function(k)
        max(abs(all_fits[[k]]$psi - ref)) > audit_par_tol, logical(1)))
    }
    gain <- if (is.finite(base$logLik)) ans$logLik - base$logLik else NA_real_
    ans$branch_ambiguous <- branch
    ans$start_discordance <- is.finite(gain) && gain > avj$effective_ll_tol
    ans$audit_ok <- TRUE
    ans$audit <- list(
      fits = all_fits, eligible = eligible, best = best,
      branch_ambiguous = branch,
      continuation_logLik_gain = gain,
      effective_ll_tol = avj$effective_ll_tol
    )
    cache_put(eta_grid[j], ans)
    update_row(j, ans, audited = TRUE)
    invisible(ans)
  }

  c_idx <- which.min(abs(eta_grid - center_eta))
  warm <- fit_single(eta_grid[c_idx], center_psi, use_cache = FALSE)
  update_row(c_idx, warm, audited = FALSE)
  warm_psi <- if (!is.null(warm$psi_hat) && all(is.finite(warm$psi_hat)))
    warm$psi_hat else center_psi
  if (c_idx < ngrid) {
    for (j in (c_idx + 1L):ngrid) {
      fc <- fit_single(eta_grid[j], warm_psi, use_cache = FALSE)
      update_row(j, fc, audited = FALSE)
      if (!is.null(fc$psi_hat) && all(is.finite(fc$psi_hat)))
        warm_psi <- fc$psi_hat
    }
  }
  warm_psi <- if (!is.null(fits[[c_idx]]$psi_hat) &&
                  all(is.finite(fits[[c_idx]]$psi_hat)))
    fits[[c_idx]]$psi_hat else center_psi
  if (c_idx > 1L) {
    for (j in (c_idx - 1L):1L) {
      fc <- fit_single(eta_grid[j], warm_psi, use_cache = FALSE)
      update_row(j, fc, audited = FALSE)
      if (!is.null(fc$psi_hat) && all(is.finite(fc$psi_hat)))
        warm_psi <- fc$psi_hat
    }
  }

  selected_p <- function() tab$p_chisq
  audit_indices <- function() {
    if (audit_spec$mode == "none") return(integer(0))
    if (audit_spec$mode == "all") return(seq_len(ngrid))
    p <- selected_p()
    acc <- tab$resolved & is.finite(p) & p >= alpha
    trans <- which(diff(as.integer(acc)) != 0L)
    unresolved <- which(!tab$resolved)
    un_neigh <- unique(c(unresolved - 1L, unresolved, unresolved + 1L))
    checkpoints <- seq.int(1L, ngrid, by = audit_spec$stride)
    required_eta <- if (length(audit_eta_values)) {
      unique(vapply(audit_eta_values, function(e) which.min(abs(eta_grid - e)), integer(1)))
    } else integer(0)
    unique(pmin(pmax(c(
      1L, ngrid, c_idx, which.min(abs(eta_grid)), required_eta,
      checkpoints, trans, trans + 1L, un_neigh,
      which(tab$reference_dominated)
    ), 1L), ngrid))
  }

  ai <- audit_indices()
  if (length(ai)) for (j in ai) audit_one(j)
  # A second, inexpensive pass catches a transition created by an audited fit.
  if (audit_spec$mode == "adaptive") {
    ai2 <- setdiff(audit_indices(), ai)
    if (length(ai2)) for (j in ai2) audit_one(j)
  }

  profile_reference_dominated <- any(tab$reference_dominated, na.rm = TRUE)
  profile_ok <- profile_reference_ok && !profile_reference_dominated
  if (!profile_ok) tab$profile_resolved[] <- FALSE

  p_sel <- selected_p()
  accept_score <- tab$resolved & is.finite(p_sel) & p_sel >= alpha
  accept_plr <- profile_ok & tab$profile_resolved &
    is.finite(tab$LR) & tab$LR <= cutoff

  margin_fit <- function(eta, warm) fit_single(eta, warm, use_cache = TRUE)
  margin_fun_score <- function(eta, warm_psi) {
    fc <- margin_fit(eta, warm_psi)
    if (!isTRUE(fc$score_usable) || !is.finite(fc$T_score)) return(NA_real_)
    cutoff - fc$T_score
  }
  margin_fun_plr <- function(eta, warm_psi) {
    fc <- margin_fit(eta, warm_psi)
    if (!isTRUE(fc$optim_ok) || !isTRUE(fc$likelihood_ok)) return(NA_real_)
    raw <- 2 * (ll_ref - fc$logLik)
    if (!is.finite(raw) || raw < -profile_lr_tol) return(NA_real_)
    cutoff - max(0, raw)
  }

  runs_to_intervals <- function(accept, resolved, refine_margin = NULL) {
    info <- grid_intervals_with_edges(eta_grid, accept)
    if (nrow(info$intervals) == 0L) {
      info$intervals_refined <- info$intervals
    } else if (isTRUE(refine) && !is.null(refine_margin)) {
      ref <- lapply(info$index_runs, function(ii)
        .bpss_refine_run(min(ii), max(ii), eta_grid, accept, resolved,
                         psi_list, refine_margin, tol = refine_tol))
      info$intervals_refined <- data.frame(
        lower = vapply(ref, function(z) unname(z["lower"]), numeric(1)),
        upper = vapply(ref, function(z) unname(z["upper"]), numeric(1)))
    } else {
      info$intervals_refined <- info$intervals
    }
    rho_refined <- transform_eta_intervals_to_rho(info$intervals_refined)
    info$length_eta <- interval_total_length(info$intervals_refined)
    info$length_rho <- interval_total_length(rho_refined)
    info$length <- info$length_rho
    info$globally_resolved <- all(resolved)
    info$length_reportable <- info$globally_resolved && !info$touches_edge
    info$reportable <- if (info$length_reportable) info$length_rho else NA_real_
    info$disconnected <- length(info$index_runs) > 1L
    info$disconnected_resolved <- info$disconnected
    info$rho_intervals_refined <- rho_refined
    info
  }

  score_info <- runs_to_intervals(accept_score, tab$resolved,
                                  margin_fun_score)
  plr_info <- if (profile_ok)
    runs_to_intervals(accept_plr, tab$profile_resolved,
                      margin_fun_plr) else NULL
  score_eta <- score_info$intervals_refined
  plr_eta <- if (!is.null(plr_info)) plr_info$intervals_refined else NULL

  list(
    table = tab, fits = if (keep_fits) fits else NULL,
    cutoff = cutoff, alpha = alpha, level = level,
    audit_mode = audit_spec$mode,
    audit_stride = audit_spec$stride,
    audit_eta_values = audit_eta_values,
    n_audited = sum(tab$audit_performed),
    profile_requested = profile_requested,
    profile_reference_ok = profile_reference_ok,
    profile_reference_audited = !is.null(profile_reference) &&
      isTRUE(profile_reference$audit_ok),
    profile_reference_dominated = profile_reference_dominated,
    profile_ok = profile_ok,
    score_eta_set = score_eta,
    score_rho_set = score_info$rho_intervals_refined,
    score_info = score_info,
    profile_eta_set = plr_eta,
    profile_rho_set = if (!is.null(plr_info))
      plr_info$rho_intervals_refined else NULL,
    profile_info = plr_info,
    n_resolved = sum(tab$resolved),
    n_unresolved = sum(!tab$resolved),
    unresolved_rho_ranges = .bpss_unresolved_ranges(tab),
    globally_resolved = all(tab$resolved),
    n_branch_ambiguous = sum(tab$branch_ambiguous, na.rm = TRUE),
    n_start_discordant = sum(tab$start_discordance, na.rm = TRUE),
    disconnected = score_info$disconnected,
    edge = score_info$touches_edge
  )
}

# Likelihood-ratio test of rho = 0 against the (audited) unrestricted fit.
bpss_lr_ignorability <- function(
    D, XS, XO, unrestricted, null_fit = NULL,
    require_audited_reference = TRUE, exclude_boundary = TRUE,
    level = 0.95, lr_tol = 1e-3,
    eps = 1e-12, rcond_tol = 1e-12,
    coef_bound = 12, maxit = 1000, optim_factr = 1e7,
    optim_pgtol = 1e-8, stationarity_tol = 1e-4,
    boundary_fraction = 0.99, ...) {
  if (is.null(null_fit))
    null_fit <- bpss_fit_null_factorized(
      D, XS, XO, eps = eps, rcond_tol = rcond_tol,
      coef_bound = coef_bound, maxit = maxit,
      optim_factr = optim_factr, optim_pgtol = optim_pgtol,
      stationarity_tol = stationarity_tol,
      boundary_fraction = boundary_fraction, check_generic = TRUE
    )

  ref_numerical_ok <- !is.null(unrestricted) &&
    isTRUE(unrestricted$optim_ok) && isTRUE(unrestricted$likelihood_ok) &&
    is.finite(unrestricted$logLik)
  ref_audit_ok <- isTRUE(unrestricted$audit_ok)
  reference_ok <- ref_numerical_ok &&
    (!isTRUE(require_audited_reference) || ref_audit_ok)
  null_ok <- isTRUE(null_fit$optim_ok) &&
             isTRUE(null_fit$likelihood_ok) &&
             is.finite(null_fit$logLik)

  LR_raw <- if (reference_ok && null_ok)
    2 * (unrestricted$logLik - null_fit$logLik) else NA_real_
  reference_dominated <- is.finite(LR_raw) && LR_raw < -lr_tol
  LR <- if (is.finite(LR_raw) && !reference_dominated)
    max(0, LR_raw) else NA_real_
  boundary_ok <- !isTRUE(exclude_boundary) ||
    (!isTRUE(unrestricted$eta_boundary_hit) &&
     !isTRUE(unrestricted$coef_boundary_hit) &&
     !isTRUE(null_fit$coef_boundary_hit))
  usable <- reference_ok && null_ok && !reference_dominated &&
            boundary_ok && is.finite(LR)
  p <- if (usable) stats::pchisq(LR, df = 1, lower.tail = FALSE) else NA_real_

  list(
    ok = usable, usable = usable, statistic = LR, statistic_raw = LR_raw,
    p_chisq = p, reject = if (usable) p < 1 - level else NA,
    reference_ok = reference_ok, reference_audited = ref_audit_ok,
    null_ok = null_ok, boundary_ok = boundary_ok,
    reference_dominated = reference_dominated,
    unrestricted = unrestricted, null = null_fit
  )
}

# Wald interval eta_hat +/- z / sqrt(J) from the unrestricted fit, with J the
# expected efficient information; the rho interval is its tanh image. It is
# computed only when the unrestricted fit passes its checks and no parameter is
# near its bound.
bpss_wald <- function(object, level = 0.95) {
  fit <- object$unrestricted
  if (is.null(fit) || !isTRUE(fit$wald_usable))
    return(list(
      eta = data.frame(lower = numeric(0), upper = numeric(0)),
      rho = data.frame(lower = numeric(0), upper = numeric(0)),
      se_eta = NA_real_, se_rho = NA_real_, level = level,
      usable = FALSE,
      reason = if (is.null(fit)) "unrestricted fit unavailable"
               else if (!isTRUE(fit$ok)) fit$status
               else if (isTRUE(fit$eta_boundary_hit))
                 "the estimate of eta is near its bound"
               else if (isTRUE(fit$coef_boundary_hit))
                 "a regression coefficient is near its bound"
               else fit$status
    ))
  eta_int <- wald_eta_interval(fit$eta_hat, fit$J, level = level)
  rho_hat <- tanh(fit$eta_hat)
  list(
    eta = eta_int,
    rho = transform_eta_intervals_to_rho(eta_int),
    se_eta = sqrt(1 / fit$J),
    se_rho = (1 - rho_hat^2) / sqrt(fit$J),
    level = level, usable = TRUE
  )
}

# Score-test p-value curve on the stored grid.
bpss_score_curve <- function(object) {
  gi <- object$grid
  if (is.null(gi))
    stop("no grid inversion stored; refit with bpss_control(do_grid = TRUE).")
  tab <- gi$table
  data.frame(eta = tab$eta, rho = tab$rho, T = tab$T, p = tab$p_chisq)
}

#' Relative profile likelihood on the inversion grid
#'
#' Returns the relative profile likelihood
#' \eqn{R_n(\eta^\star) = \exp\{\ell_n^p(\eta^\star) - \sup_\eta \ell_n^p(\eta)\}}{R(eta*) = exp{lp(eta*) - sup lp(eta)}}
#' at the grid values of a fitted \code{"bpss"} object, where
#' \eqn{\ell_n^p}{lp} is the profile log-likelihood and the multi-start
#' unrestricted fit is the reference, together with the sets on which
#' \eqn{R_n}{R} is at least each of \code{levels}.
#'
#' @param object A fitted object of class \code{"bpss"}.
#' @param levels Relative-likelihood levels defining the reported sets.
#'
#' @return A list with \code{curve}, a data frame with columns \code{eta},
#'   \code{rho}, \code{LR} (\eqn{-2\log R_n}{-2 log R}) and \code{R}
#'   (\eqn{R_n}{R}), \code{levels}, and \code{rho_sets}, the sets on the
#'   correlation scale at grid resolution.
#' @seealso \code{\link{bpss}}, \code{\link{plot.bpss}}
#' @export
bpss_profile <- function(object, levels = c(1/8, 1/32)) {
  if (!inherits(object, "bpss")) stop("object must be of class 'bpss'.")
  gi <- object$grid
  if (is.null(gi) || !isTRUE(gi$profile_ok))
    stop("the relative profile likelihood requires a grid inversion with an ",
         "audited unrestricted fit that no constrained fit exceeds.")
  tab <- gi$table
  cuts <- -2 * log(levels)
  curve <- data.frame(eta = tab$eta, rho = tab$rho, LR = tab$LR,
                      R = exp(-tab$LR / 2))
  sets <- lapply(cuts, function(cc) {
    acc <- is.finite(tab$LR) & (tab$LR <= cc)
    transform_eta_intervals_to_rho(grid_intervals_with_edges(tab$eta, acc)$intervals)
  })
  names(sets) <- paste0("level_", format(levels))
  list(curve = curve, levels = levels, cuts = cuts, rho_sets = sets)
}

#' Test of ignorable selection
#'
#' Efficient-score test of \eqn{\rho = 0}{rho = 0}. At \eqn{\eta = 0}{eta = 0}
#' the constrained likelihood separates into a selection probit on all units and
#' an outcome probit on the selected units. The efficient score is evaluated at
#' that fit, with the expected efficient information
#' \eqn{\hat J_n(0)}{J(0)} as its variance. The likelihood-ratio statistic
#' against the multi-start unrestricted fit is reported for comparison.
#'
#' @param object A fitted object of class \code{"bpss"}.
#' @param require_audited_lr Logical; compute the likelihood-ratio statistic
#'   only when the unrestricted reference fit passed the multi-start search.
#'
#' @return A list with components \code{ok}, \code{status}, \code{U},
#'   \code{J}, \code{I_etaeta}, \code{J_over_I_etaeta}, \code{Delta},
#'   \code{T}, \code{p_chisq}, \code{LR}, \code{p_lr}, \code{lr_usable} and
#'   the null fit.
#' @export
test_ignorability <- function(object, require_audited_lr = TRUE) {
  if (!inherits(object, "bpss")) stop("object must be of class 'bpss'.")
  XS <- object$XS; XO <- object$XO; D <- object$D
  ctrl <- object$control

  fc <- NULL
  if (!is.null(object$grid)) {
    j0 <- which(abs(object$grid$table$eta) < 1e-10)
    if (length(j0) == 1L && !is.null(object$grid$fits))
      fc <- object$grid$fits[[j0]]
  }
  if (is.null(fc) || !isTRUE(fc$score_usable)) {
    fc <- bpss_fit_null_factorized(
      D, XS, XO, eps = ctrl$eps,
      rcond_tol = ctrl$rcond_tol,
      coef_bound = ctrl$coef_bound, maxit = ctrl$maxit,
      optim_factr = ctrl$optim_factr,
      optim_pgtol = ctrl$optim_pgtol,
      stationarity_tol = ctrl$stationarity_tol,
      boundary_fraction = ctrl$boundary_fraction,
      check_generic = TRUE,
      generic_ll_tol = ctrl$audit_ll_tol,
      generic_par_tol = ctrl$audit_par_tol
    )
  }
  if (!isTRUE(fc$score_usable))
    return(list(
      ok = FALSE, status = fc$status, eta0 = 0,
      U = NA_real_, J = NA_real_, I_etaeta = NA_real_,
      J_over_I_etaeta = NA_real_, Delta = NA_real_,
      T = NA_real_, p_chisq = NA_real_, LR = NA_real_, p_lr = NA_real_,
      lr_usable = FALSE, null_fit = fc
    ))

  out <- list(
    ok = TRUE, status = "ok", eta0 = 0, U = fc$U, J = fc$J,
    I_etaeta = fc$I_etaeta,
    J_over_I_etaeta = fc$J / fc$I_etaeta,
    Delta = fc$Delta, T = fc$T_score,
    p_chisq = stats::pchisq(fc$T_score, df = 1, lower.tail = FALSE),
    boundary_hit = fc$boundary_hit,
    coef_boundary_hit = fc$coef_boundary_hit,
    factorized_null = identical(fc$fit_role, "null_factorized"),
    generic_null_agrees = fc$generic_agrees %||% NA,
    null_fit = fc
  )
  lr <- bpss_lr_ignorability(
    D, XS, XO, unrestricted = object$unrestricted,
    null_fit = fc,
    require_audited_reference = require_audited_lr,
    exclude_boundary = TRUE, level = object$level,
    lr_tol = ctrl$profile_lr_tol,
    eps = ctrl$eps, rcond_tol = ctrl$rcond_tol,
    coef_bound = ctrl$coef_bound, maxit = ctrl$maxit,
    optim_factr = ctrl$optim_factr,
    optim_pgtol = ctrl$optim_pgtol,
    stationarity_tol = ctrl$stationarity_tol,
    boundary_fraction = ctrl$boundary_fraction
  )
  out$lr <- lr
  out$LR <- lr$statistic
  out$p_lr <- lr$p_chisq
  out$lr_usable <- lr$usable
  out
}

#' Efficient-score test at specified correlations
#'
#' Computes
#' \eqn{T_n(\eta^\star) = \hat U_n(\eta^\star)^2 / \hat J_n(\eta^\star)}{T(eta*) = U(eta*)^2 / J(eta*)}
#' at each hypothesized correlation \eqn{\rho^\star = \tanh(\eta^\star)}{rho* = tanh(eta*)},
#' where \eqn{\hat U_n}{U} is the efficient score and \eqn{\hat J_n}{J} the
#' expected efficient information, both evaluated at the regression
#' coefficients that maximize the likelihood with \eqn{\eta}{eta} fixed at
#' \eqn{\eta^\star}{eta*}. The p-value is computed from the
#' \eqn{\chi^2_1}{chi-squared(1)} distribution.
#'
#' The constrained fit is started from the unrestricted estimate, from separate
#' selection and selected-outcome probit fits and, when it differs, from the
#' single-start unrestricted fit; the fit with the largest log-likelihood is
#' used. The statistic is not computed if this fit fails its convergence,
#' stationarity, cell-probability or nuisance-information checks, or if a
#' regression coefficient is near its bound.
#'
#' @param object A fitted object of class \code{"bpss"}.
#' @param rho Hypothesized correlations in \eqn{(-1, 1)}{(-1, 1)}.
#'
#' @return A data frame with one row per value of \code{rho}, giving
#'   \code{T}, \code{p_value}, the efficient score \code{U}, \code{J},
#'   \code{I_etaeta}, \code{J_over_I_etaeta}, \code{Delta}, the constrained
#'   log-likelihood, whether the statistic was \code{computed}, and the fit
#'   status.
#' @seealso \code{\link{test_ignorability}}, \code{\link{confint.bpss}}
#' @export
bpss_score_test <- function(object, rho) {
  if (!inherits(object, "bpss")) stop("object must be of class 'bpss'.")
  rho <- as.numeric(rho)
  if (!length(rho) || any(!is.finite(rho)) || any(abs(rho) >= 1))
    stop("rho must contain values in (-1, 1).")
  ctrl <- object$control
  usable_psi <- function(f) !is.null(f) && isTRUE(f$likelihood_ok) &&
    !is.null(f$psi) && all(is.finite(f$psi))
  warm <- if (usable_psi(object$unrestricted)) object$unrestricted$psi else NULL
  extras <- if (usable_psi(object$unrestricted_default))
    list(object$unrestricted_default$psi) else NULL
  rows <- lapply(rho, function(r) {
    fc <- bpss_fit_constrained_audited(
      atanh(r), object$D, object$XS, object$XO,
      start_psi = warm, extra_starts = extras,
      eps = ctrl$eps, rcond_tol = ctrl$rcond_tol,
      coef_bound = ctrl$coef_bound, maxit = ctrl$maxit,
      optim_factr = ctrl$optim_factr, optim_pgtol = ctrl$optim_pgtol,
      stationarity_tol = ctrl$stationarity_tol,
      boundary_fraction = ctrl$boundary_fraction,
      ll_tol = ctrl$audit_ll_tol, par_tol = ctrl$audit_par_tol
    )
    computed <- isTRUE(fc$score_usable) && isTRUE(fc$audit_ok) &&
      is.finite(fc$T_score)
    data.frame(
      rho = r, eta = atanh(r),
      T = if (computed) fc$T_score else NA_real_,
      p_value = if (computed)
        stats::pchisq(fc$T_score, df = 1, lower.tail = FALSE) else NA_real_,
      U = fc$U, J = fc$J, I_etaeta = fc$I_etaeta,
      J_over_I_etaeta = fc$J / fc$I_etaeta, Delta = fc$Delta,
      logLik = fc$logLik, computed = computed, status = fc$status,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
