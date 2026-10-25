# S3 methods, printing, plotting, the confidence-set class and the comparison
# table across specifications.

#' Methods for fitted bpss models
#'
#' Printing, summary and extractor methods for objects returned by
#' \code{\link{bpss}}. Coefficients are displayed in the order outcome
#' (\eqn{\beta}{beta}), selection (\eqn{\gamma}{gamma}),
#' \eqn{\eta = \mathrm{atanh}(\rho)}{eta = atanh(rho)}. Standard errors use the
#' inverse expected information at the unrestricted estimate. z statistics and
#' p-values for the regression coefficients are reported only when the Wald
#' interval is computed (see \code{\link{bpss}}); none is reported for
#' \eqn{\eta}{eta}, for which the efficient-score results should be used.
#'
#' \code{summary()} reports the efficient-score confidence set with the grid on
#' which it was computed and any grid values at which the statistic was not
#' computed. For each connected component of the set it also reports whether
#' the component contains the unrestricted estimate or reaches an end of the
#' grid, the largest relative profile likelihood over its evaluated grid
#' values, and the ranges of \eqn{\hat J_n}{J},
#' \eqn{\hat J_n / \hat I_{\eta\eta,n}}{J / I_etaeta} and
#' \eqn{\hat\Delta_n}{Delta} there, as recommended in the paper.
#'
#' @param x,object A fitted object of class \code{"bpss"} (or its summary,
#'   confidence set or decomposition, as appropriate).
#' @param type For \code{coef()}: which coefficients to return.
#' @param digits Number of significant digits to print.
#' @param ... Further arguments, ignored.
#' @name bpss-methods
NULL

# Display order (beta, gamma, eta) and labels from the internal order.
.bpss_display <- function(object) {
  dd <- object$dims
  perm <- c(dd$idx_beta, dd$idx_gamma, dd$idx_eta)
  nm_int <- object$coef_names
  nm <- nm_int
  nm[dd$idx_gamma] <- sub("^gamma:", "selection:", nm_int[dd$idx_gamma])
  nm[dd$idx_beta]  <- sub("^beta:",  "outcome:",   nm_int[dd$idx_beta])
  list(perm = perm, names = nm[perm])
}

#' @rdname bpss-methods
#' @export
coef.bpss <- function(object, type = c("all", "outcome", "selection", "eta", "rho"), ...) {
  type <- match.arg(type)
  fit <- object$unrestricted
  numerical_ok <- !is.null(fit) && isTRUE(fit$optim_ok) &&
                  isTRUE(fit$likelihood_ok) &&
                  all(is.finite(fit$theta))
  if (!numerical_ok) return(NA_real_)
  if (type == "rho") return(c(rho = unname(fit$rho_hat)))
  if (type == "eta") return(c(eta = unname(fit$eta_hat)))
  dsp <- .bpss_display(object)
  b <- as.numeric(fit$theta)[dsp$perm]; names(b) <- dsp$names
  dd <- object$dims
  if (type == "outcome")   return(b[seq_along(dd$idx_beta)])
  if (type == "selection") return(b[length(dd$idx_beta) + seq_along(dd$idx_gamma)])
  b
}

#' @rdname bpss-methods
#' @export
vcov.bpss <- function(object, ...) {
  fit <- object$unrestricted
  dsp <- .bpss_display(object)
  p <- length(dsp$perm)
  empty <- matrix(NA_real_, p, p,
                  dimnames = list(dsp$names, dsp$names))
  if (is.null(fit) || !isTRUE(fit$info_ok) || is.null(fit$info))
    return(empty)
  I <- fit$info$I
  V <- tryCatch(chol2inv(chol(I)), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(V))) return(empty)
  V <- V[dsp$perm, dsp$perm, drop = FALSE]
  dimnames(V) <- list(dsp$names, dsp$names)
  V
}

#' @rdname bpss-methods
#' @export
logLik.bpss <- function(object, ...) {
  fit <- object$unrestricted
  val <- if (!is.null(fit)) fit$logLik else NA_real_
  structure(val, df = object$dims$ptheta, nobs = nrow(object$D), class = "logLik")
}

#' @rdname bpss-methods
#' @export
nobs.bpss <- function(object, ...) nrow(object$D)

# Confidence-set class: a data frame of intervals (lower, upper) with the
# level, method and diagnostic notes as attributes.
new_bpss_set <- function(intervals, parm, level, method, refined = NA,
                         disconnected = FALSE, edge = FALSE,
                         unresolved = 0L, usable = TRUE,
                         reason = NA_character_, grid_note = NA_character_,
                         unresolved_note = NA_character_) {
  if (is.null(intervals)) intervals <- data.frame(lower = numeric(0), upper = numeric(0))
  intervals <- as.data.frame(intervals)
  attr(intervals, "parm") <- parm
  attr(intervals, "level") <- level
  attr(intervals, "method") <- method
  attr(intervals, "refined") <- refined
  attr(intervals, "disconnected") <- isTRUE(disconnected)
  attr(intervals, "edge") <- isTRUE(edge)
  attr(intervals, "unresolved") <- as.integer(unresolved)
  attr(intervals, "usable") <- isTRUE(usable)
  attr(intervals, "reason") <- reason
  attr(intervals, "grid_note") <- grid_note
  attr(intervals, "unresolved_note") <- unresolved_note
  class(intervals) <- c("bpss_set", "data.frame")
  intervals
}

# Grid values at which the score statistic was not computed, as text, with
# the reasons; NULL when the statistic was computed everywhere.
.bpss_unresolved_note <- function(tab, digits = 3) {
  ii <- which(!tab$resolved)
  if (!length(ii)) return(NULL)
  runs <- split(ii, cumsum(c(TRUE, diff(ii) != 1L)))
  where <- vapply(runs, function(j)
    if (length(j) == 1L) sprintf("%.*f", digits, tab$rho[j])
    else sprintf("%.*f to %.*f", digits, min(tab$rho[j]), digits, max(tab$rho[j])),
    character(1))
  reason <- vapply(ii, function(j) {
    e <- tab$information_error[j]
    if (!is.na(e)) return(sub("[.]$", "", e))
    if (isTRUE(tab$coef_boundary_hit[j]))
      return("a regression coefficient is near its bound")
    st <- tab$status[j]
    if (is.na(st)) st <- ""
    switch(st,
           constrained_nonconvergence = "the constrained fit did not converge",
           constrained_nonstationary = "the constrained fit failed the stationarity check",
           bad_likelihood = "the log-likelihood is not finite",
           constrained_error = , constrained_optim_error = "the constrained fit failed",
           bad_information_or_score = "the information or score could not be computed",
           "the constrained fit failed its checks")
  }, character(1))
  sprintf("statistic not computed at rho = %s (%s); whether %s to the set is not determined",
          .bpss_list_text(where), paste(unique(reason), collapse = "; "),
          if (length(ii) == 1L) "this value belongs" else "these values belong")
}

.bpss_list_text <- function(x) {
  if (length(x) <= 1L) return(x)
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
}

.bpss_grid_note <- function(tab, refined) {
  sprintf("computed on %d grid values of rho in [%.3f, %.3f]%s", nrow(tab),
          min(tab$rho), max(tab$rho),
          if (isTRUE(refined)) ", with endpoints refined where the decision changes"
          else ", with endpoints at grid resolution")
}

# Three significant digits, keeping trailing zeros; scientific when small.
.fmt3 <- function(x) {
  if (length(x) != 1L || !is.finite(x)) return("NA")
  if (x != 0 && (abs(x) < 1e-3 || abs(x) >= 1e5))
    return(formatC(x, digits = 2, format = "e"))
  sub("[.]$", "", trimws(formatC(signif(x, 3), digits = 3, format = "fg", flag = "#")))
}

#' @rdname bpss-methods
#' @export
print.bpss_set <- function(x, digits = 4, ...) {
  lev <- attr(x, "level"); method <- attr(x, "method"); parm <- attr(x, "parm")
  label <- if (identical(method, "wald")) "Wald" else "efficient-score"
  cat(sprintf("%.0f%% %s confidence set on the %s scale\n", 100 * lev, label, parm))
  if (!isTRUE(attr(x, "usable"))) {
    why <- attr(x, "reason")
    cat("  not computed", if (is.character(why) && length(why) &&
                              !is.na(why) && nzchar(why))
      paste0(": ", why) else "", "\n", sep = "")
    return(invisible(x))
  }
  y <- x; class(y) <- "data.frame"
  if (nrow(y) == 0L) cat("  empty on the evaluated grid\n")
  else print(y, row.names = FALSE, digits = digits)
  notes <- character(0)
  gn <- attr(x, "grid_note")
  if (identical(method, "score") && is.character(gn) && !is.na(gn)) notes <- c(notes, gn)
  if (isTRUE(attr(x, "disconnected"))) notes <- c(notes, "disconnected set")
  if (isTRUE(attr(x, "edge")))
    notes <- c(notes, "a component reaches an end of the grid and is reported only up to it")
  un <- attr(x, "unresolved_note")
  if (is.character(un) && !is.na(un)) notes <- c(notes, un)
  if (length(notes)) cat(paste0("  note: ", notes, "\n"), sep = "")
  invisible(x)
}

#' Confidence sets for the correlation
#'
#' Efficient-score confidence set (the default) or Wald interval for
#' \eqn{\rho}{rho} or \eqn{\eta = \mathrm{atanh}(\rho)}{eta = atanh(rho)}.
#'
#' At the level used in the fit, the score set is the stored set with refined
#' endpoints. At another level it is recomputed from the stored grid values,
#' with endpoints at grid resolution. The printed notes give the grid, any grid
#' values at which the statistic was not computed, and whether a component
#' reaches an end of the grid.
#'
#' @param object A fitted object of class \code{"bpss"}.
#' @param parm \code{"rho"} or \code{"eta"}.
#' @param level Confidence level; defaults to the level used in the fit.
#' @param method \code{"score"} or \code{"wald"}.
#' @param ... Ignored.
#' @return An object of class \code{"bpss_set"}: a data frame with columns
#'   \code{lower} and \code{upper}, one row per connected component, with the
#'   level, method and diagnostic flags as attributes.
#' @export
confint.bpss <- function(object, parm = c("rho", "eta"), level = NULL,
                         method = c("score", "wald"), ...) {
  parm <- match.arg(parm); method <- match.arg(method)
  if (is.null(level)) level <- object$level
  if (method == "wald") {
    w <- if (isTRUE(all.equal(level, object$level))) object$wald else bpss_wald(object, level = level)
    set <- if (parm == "rho") w$rho else w$eta
    return(new_bpss_set(
      set, parm = parm, level = level, method = "wald", refined = NA,
      disconnected = FALSE, edge = FALSE, unresolved = 0L,
      usable = isTRUE(w$usable), reason = w$reason %||% NA_character_
    ))
  }
  gi <- object$grid
  if (is.null(gi)) stop("no grid inversion stored; refit with control = bpss_control(do_grid = TRUE).")
  alpha <- 1 - level
  refined <- isTRUE(all.equal(level, object$level))
  if (refined) {
    set <- if (parm == "rho") gi$score_rho_set else gi$score_eta_set
    disc <- isTRUE(gi$disconnected); edge <- isTRUE(gi$edge); unres <- gi$n_unresolved
  } else {
    curve <- bpss_score_curve(object)
    acc <- gi$table$resolved & is.finite(curve$p) & (curve$p >= alpha)
    info <- grid_intervals_with_edges(gi$table$eta, acc)
    eta_int <- info$intervals
    set <- if (parm == "rho") transform_eta_intervals_to_rho(eta_int) else eta_int
    disc <- length(info$index_runs) > 1L
    edge <- isTRUE(info$touches_edge)
    unres <- sum(!(gi$table$resolved & is.finite(curve$p)))
  }
  score_usable <- gi$n_resolved > 0L
  un <- .bpss_unresolved_note(gi$table)
  new_bpss_set(
    set, parm = parm, level = level, method = "score", refined = refined,
    disconnected = disc, edge = edge, unresolved = unres,
    usable = score_usable,
    reason = if (score_usable) NA_character_
             else "the statistic was not computed at any grid value",
    grid_note = .bpss_grid_note(gi$table, refined && isTRUE(object$control$refine_grid)),
    unresolved_note = if (is.null(un)) NA_character_ else un
  )
}

#' Fitted cell probabilities
#'
#' @param object A fitted object of class \code{"bpss"}.
#' @param newdata Optional data frame with the covariates of both equations;
#'   by default the data used in the fit.
#' @param type \code{"cell"} for the probabilities of the three observed cells
#'   (S = 0; S = 1, Y = 0; S = 1, Y = 1), \code{"selection"} for
#'   \eqn{\mathrm{pr}(S = 1)}{pr(S = 1)}, or \code{"conditional"} for
#'   \eqn{\mathrm{pr}(Y = 1 \mid S = 1)}{pr(Y = 1 | S = 1)}.
#' @param ... Ignored.
#' @return A matrix (\code{type = "cell"}) or a numeric vector, evaluated at
#'   the unrestricted estimate.
#' @export
predict.bpss <- function(object, newdata = NULL,
                         type = c("cell", "selection", "conditional"), ...) {
  type <- match.arg(type)
  fit <- object$unrestricted
  numerical_ok <- !is.null(fit) && isTRUE(fit$optim_ok) &&
                  isTRUE(fit$likelihood_ok) &&
                  all(is.finite(fit$theta))
  if (!numerical_ok) stop("no numerically valid unrestricted fit.")
  if (is.null(newdata)) {
    XS <- object$XS; XO <- object$XO
  } else {
    XS <- stats::model.matrix(stats::delete.response(object$terms$selection), newdata,
                              xlev = object$xlevels$selection)
    XO <- stats::model.matrix(stats::delete.response(object$terms$outcome), newdata,
                              xlev = object$xlevels$outcome)
  }
  prob <- bpss_cell_prob(fit$theta, XS, XO)
  if (type == "cell") return(prob)
  if (type == "selection") return(prob[, "p10"] + prob[, "p11"])          # P(S = 1)
  prob[, "p11"] / pmax(prob[, "p10"] + prob[, "p11"], .Machine$double.eps) # P(Y = 1 | S = 1)
}

#' @rdname predict.bpss
#' @export
fitted.bpss <- function(object, ...) predict(object, type = "cell")

#' @rdname bpss-methods
#' @export
summary.bpss <- function(object, ...) {
  fit <- object$unrestricted
  dsp <- .bpss_display(object)
  estimable <- !is.null(fit) &&
               isTRUE(fit$likelihood_ok) &&
               all(is.finite(fit$theta))
  est <- if (estimable)
    as.numeric(fit$theta)[dsp$perm]
    else rep(NA_real_, length(dsp$perm))
  V <- vcov(object)
  se <- if (all(dim(V) > 0) && isTRUE(fit$info_ok))
    sqrt(diag(V)) else rep(NA_real_, length(est))
  wald_valid <- isTRUE(fit$wald_usable)
  z <- if (wald_valid) est / se else rep(NA_real_, length(est))
  p <- if (wald_valid)
    2 * stats::pnorm(abs(z), lower.tail = FALSE)
  else rep(NA_real_, length(est))
  # No Wald test for eta: inference on rho uses the efficient score.
  is_eta <- dsp$names == "eta"
  z[is_eta] <- NA_real_
  p[is_eta] <- NA_real_
  coef_tab <- cbind(
    Estimate = est, `Std. Error` = se,
    `z value` = z, `Pr(>|z|)` = p
  )
  rownames(coef_tab) <- dsp$names

  rho_hat <- if (!is.null(fit)) fit$rho_hat else NA_real_
  se_rho <- if (!is.null(fit) && is.finite(fit$J) && fit$J > 0)
    (1 - rho_hat^2) / sqrt(fit$J) else NA_real_
  out <- list(
    call = object$call, coefficients = coef_tab,
    rho_hat = rho_hat, se_rho = se_rho,
    eta_hat = if (!is.null(fit)) fit$eta_hat else NA_real_,
    J = if (!is.null(fit)) fit$J else NA_real_,
    I_etaeta = if (!is.null(fit)) fit$I_etaeta else NA_real_,
    Delta = if (!is.null(fit)) fit$Delta else NA_real_,
    status = if (!is.null(fit)) fit$status else NA_character_,
    optim_ok = if (!is.null(fit)) fit$optim_ok else FALSE,
    stationary = if (!is.null(fit)) fit$stationary else FALSE,
    projected_gradient_norm =
      if (!is.null(fit)) fit$projected_gradient_norm else NA_real_,
    boundary_hit = if (!is.null(fit)) fit$boundary_hit else NA,
    eta_boundary_hit =
      if (!is.null(fit)) fit$eta_boundary_hit else NA,
    coef_boundary_hit =
      if (!is.null(fit)) fit$coef_boundary_hit else NA,
    wald_usable = if (!is.null(fit)) fit$wald_usable else FALSE,
    audit_ok = if (!is.null(fit)) fit$audit_ok %||% NA else NA,
    default_discordance =
      if (!is.null(fit)) fit$default_discordance %||% NA else NA,
    logLik = if (!is.null(fit)) fit$logLik else NA_real_,
    nobs = nrow(object$D),
    n_selected = sum(object$S == 1L),
    excluded = object$excluded, level = object$level,
    wald = object$wald, grid = object$grid,
    ignorability = object$ignorability,
    components = .bpss_components(object),
    grid_refined = isTRUE(object$control$refine_grid),
    dims = object$dims
  )
  class(out) <- "summary.bpss"
  out
}

.fmt_set <- function(set, digits = 3) {
  if (is.null(set) || nrow(set) == 0L) return("empty")
  paste(apply(set, 1L, function(r) sprintf("[%.*f, %.*f]", digits, r[1], digits, r[2])), collapse = " U ")
}

# Diagnostics recommended in the paper for each connected component of the score set:
# refined endpoints (rho scale), whether the component contains the
# unrestricted estimate or reaches a grid end, the largest relative profile
# likelihood over its grid values and the ranges of J, J / I_etaeta and Delta.
.bpss_components <- function(object) {
  gi <- object$grid
  if (is.null(gi) || is.null(gi$score_info)) return(NULL)
  info <- gi$score_info
  runs <- info$index_runs
  if (!length(runs)) return(NULL)
  tab <- gi$table
  ngrid <- nrow(tab)
  fit <- object$unrestricted
  rho_hat <- if (!is.null(fit) && isTRUE(fit$optim_ok) &&
                 isTRUE(fit$likelihood_ok) && is.finite(fit$eta_hat))
    tanh(fit$eta_hat) else NA_real_
  R <- if (isTRUE(gi$profile_ok)) exp(-tab$LR / 2) else rep(NA_real_, ngrid)
  rng <- function(v) {
    v <- v[is.finite(v)]
    if (length(v)) range(v) else c(NA_real_, NA_real_)
  }
  out <- do.call(rbind, lapply(seq_along(runs), function(k) {
    ii <- runs[[k]]
    lo <- tanh(info$intervals_refined$lower[k])
    hi <- tanh(info$intervals_refined$upper[k])
    J <- rng(tab$J[ii]); ratio <- rng(tab$J[ii] / tab$I_etaeta[ii])
    Delta <- rng(tab$Delta[ii])
    data.frame(
      lower = lo, upper = hi,
      contains_estimate = is.finite(rho_hat) && rho_hat >= lo && rho_hat <= hi,
      reaches_boundary = min(ii) == 1L || max(ii) == ngrid,
      max_R = if (any(is.finite(R[ii]))) max(R[ii], na.rm = TRUE) else NA_real_,
      J_min = J[1L], J_max = J[2L],
      J_ratio_min = ratio[1L], J_ratio_max = ratio[2L],
      Delta_min = Delta[1L], Delta_max = Delta[2L]
    )
  }))
  rownames(out) <- NULL
  out
}

#' @rdname bpss-methods
#' @export
print.summary.bpss <- function(x, digits = 4, ...) {
  cat("Bivariate probit sample selection model\n")
  cat("Call: "); print(x$call)
  cat(sprintf("\nObservations: %d (selected: %d)\n",
              x$nobs, x$n_selected))
  cat("\nCoefficients (outcome, selection, eta = atanh(rho)):\n")
  stats::printCoefmat(
    x$coefficients, digits = digits, signif.stars = FALSE,
    has.Pvalue = TRUE, na.print = ""
  )
  if (isTRUE(x$wald_usable))
    cat("No z statistic is given for eta; use the efficient-score results for rho.\n")
  else
    cat(sprintf("No z statistics are given: %s.\n",
                x$wald$reason %||% "the unrestricted fit failed its checks"))
  cat(sprintf(
    "\nrho_hat = %.*f   eta_hat = %.*f   largest projected-gradient component = %.*g\n",
    digits, x$rho_hat, digits, x$eta_hat,
    max(digits, 4L), x$projected_gradient_norm
  ))
  if (is.finite(x$J))
    cat(sprintf(
      "At the estimate: J_n = %.*f; J_n / I_etaeta,n = %.*g; Delta_n = %.*g\n",
      digits, x$J, max(digits, 4L), x$J / x$I_etaeta,
      max(digits, 4L), x$Delta
    ))
  else
    cat("At the estimate: the efficient information was not computed.\n")
  if (isTRUE(x$eta_boundary_hit))
    cat("Note: the estimate of eta is near its bound.\n")
  if (isTRUE(x$coef_boundary_hit))
    cat("Note: a regression coefficient is near its bound.\n")
  if (isTRUE(x$audit_ok))
    cat("Unrestricted fit: multi-start search.\n")
  if (isTRUE(x$default_discordance))
    cat("Note: the single-start and multi-start fits differ materially under the audit rule.\n")

  if (!is.null(x$grid)) {
    tab <- x$grid$table
    cat(sprintf("\n%.0f%% efficient-score set (rho): %s\n",
                100 * x$level, .fmt_set(x$grid$score_rho_set, digits = 3)))
    cat(sprintf("  %s.\n", .bpss_grid_note(tab, x$grid_refined)))
    un <- .bpss_unresolved_note(tab)
    if (!is.null(un))
      cat(sprintf("  %s%s.\n", toupper(substr(un, 1L, 1L)), substring(un, 2L)))
    if (isTRUE(x$grid$edge))
      cat("  A component reaches an end of the grid and is reported only up to it.\n")
    if (x$grid$n_branch_ambiguous > 0L)
      cat(sprintf("  %d grid value(s) with distinct nuisance fits of equal likelihood.\n",
                  x$grid$n_branch_ambiguous))
    if (isTRUE(x$grid$profile_reference_dominated))
      cat("  A constrained fit exceeds the unrestricted reference; the profile likelihood is not reported.\n")
    cm <- x$components
    if (!is.null(cm) && nrow(cm)) {
      cat("\nComponents of the score set (over their evaluated grid values):\n")
      for (k in seq_len(nrow(cm))) {
        flags <- c(if (isTRUE(cm$contains_estimate[k])) "contains the estimate",
                   if (isTRUE(cm$reaches_boundary[k])) "reaches an end of the grid")
        cat(sprintf("  [%.3f, %.3f]%s\n", cm$lower[k], cm$upper[k],
                    if (length(flags)) paste0(": ", paste(flags, collapse = "; ")) else ""))
        cat(sprintf("    largest relative profile likelihood %s\n",
                    if (is.finite(cm$max_R[k])) .fmt3(cm$max_R[k]) else "not available"))
        cat(sprintf("    J_n %s to %s; J_n / I_etaeta,n %s to %s; Delta_n %s to %s\n",
                    .fmt3(cm$J_min[k]), .fmt3(cm$J_max[k]),
                    .fmt3(cm$J_ratio_min[k]), .fmt3(cm$J_ratio_max[k]),
                    .fmt3(cm$Delta_min[k]), .fmt3(cm$Delta_max[k])))
      }
    }
  }

  if (!is.null(x$wald)) {
    if (isTRUE(x$wald$usable))
      cat(sprintf("\n%.0f%% Wald interval (rho), for comparison: %s\n",
                  100 * x$level,
                  .fmt_set(x$wald$rho, digits = 3)))
    else
      cat(sprintf("\nWald interval not computed: %s.\n",
                  x$wald$reason %||% "the unrestricted fit failed its checks"))
  }

  ig <- x$ignorability
  if (!is.null(ig) && isTRUE(ig$ok)) {
    cat(sprintf(
      "\nTest of ignorable selection (rho = 0): score T = %.3f, p = %.4f",
      ig$T, ig$p_chisq
    ))
    if (isTRUE(ig$lr_usable)) {
      cat(sprintf("; LR = %.3f, p = %.4f", ig$LR, ig$p_lr))
      if (isTRUE(ig$lr_reference_at_bound))
        cat(sprintf(" (unrestricted estimate near the bound |rho| = %.3f)", ig$rho_bound))
    } else {
      cat(sprintf("; LR not computed (%s)", ig$lr_reason %||% "no valid reference fit"))
    }
    cat(sprintf(
      "\n  J_n(0) = %.*f; J_n(0) / I_etaeta,n(0) = %.*g; Delta_n(0) = %.*g\n",
      digits, ig$J, max(digits, 4L),
      ig$J_over_I_etaeta, max(digits, 4L), ig$Delta
    ))
    if (identical(ig$generic_null_agrees, FALSE))
      cat("  Note: the generic constrained fit at rho = 0 did not reproduce the factorized fit.\n")
  } else if (!is.null(ig)) {
    cat(sprintf("\nTest of ignorable selection not computed: %s.\n",
                ig$status %||% "the constrained fit at rho = 0 failed"))
  }
  if (length(x$excluded))
    cat(sprintf("\nSelection-only covariate(s): %s\n",
                paste(x$excluded, collapse = ", ")))
  invisible(x)
}

#' @rdname bpss-methods
#' @export
print.bpss <- function(x, digits = 4, ...) {
  fit <- x$unrestricted
  cat("Bivariate probit sample selection model\n")
  cat("Call: "); print(x$call)
  numerical_fit <- !is.null(fit) &&
                   isTRUE(fit$likelihood_ok) &&
                   !is.null(fit$theta) &&
                   all(is.finite(fit$theta))
  if (!numerical_fit) {
    cat("\n** unrestricted fit unavailable (status: ",
        if (!is.null(fit)) fit$status else "NA", ") **\n", sep = "")
    return(invisible(x))
  }
  if (!isTRUE(fit$stationary))
    cat("\n** unrestricted numerical fit failed the stationarity check **\n")
  cat(sprintf(
    "\nrho_hat = %.*f   eta_hat = %.*f   J = %.*f   Delta = %.*g\n",
    digits, fit$rho_hat, digits, fit$eta_hat,
    digits, fit$J, max(digits, 4L), fit$Delta
  ))
  if (!is.null(x$grid))
    cat(sprintf(
      "%.0f%% efficient-score set (rho): %s\n",
      100 * x$level,
      .fmt_set(x$grid$score_rho_set)
    ))
  if (!is.null(x$ignorability) && isTRUE(x$ignorability$ok)) {
    cat(sprintf("Test of rho = 0: score p = %.4f\n",
                x$ignorability$p_chisq))
    if (isTRUE(x$ignorability$lr_usable))
      cat(sprintf("Test of rho = 0: LR p = %.4f\n",
                  x$ignorability$p_lr))
  }
  cat("Use summary() for coefficients and diagnostics.\n")
  invisible(x)
}

#' @rdname bpss-methods
#' @export
print.bpss_decomp <- function(x, digits = 4, ...) {
  where <- switch(x$at %||% "supplied",
                  mle = "the unrestricted estimate",
                  null = "the constrained fit at rho = 0",
                  "the supplied parameter values")
  cat(sprintf("Local exclusion decomposition for '%s', evaluated at %s\n",
              x$exclusion, where))
  cat(sprintf("  gamma_V = %.4g (%s coordinates)\n", x$gamma_V,
              if (isTRUE(x$residualize)) "residualized" else "unresidualized"))
  cat(sprintf("  K0 = %.5g   K1 = %.5g   K01 = %.5g\n", x$K0, x$K1, x$K01))
  cat(sprintf("  baseline K0 = %.5g; exclusion term gamma_V^2 K1 = %.5g; cross term 2 gamma_V K01 = %.5g\n",
              x$K0, x$exclusion_term, x$cross_term))
  if (is.finite(x$kappa)) {
    sgn <- if (x$kappa > 0) "positive local cross term"
           else if (x$kappa < 0) "negative local cross term"
           else "zero local cross term"
    cat(sprintf("  kappa = %.4f (%s)   omega = %.4g\n", x$kappa, sgn, x$omega))
  } else {
    cat(sprintf("  kappa and omega not reported: %s\n",
                x$kappa_omega_undefined_reason))
  }
  ratio <- x$B / x$J_over_n
  cat(sprintf("  A = %.5g   B = %.5g   J/n = %.5g   B/(J/n) = %s\n",
              x$A, x$B, x$J_over_n, .fmt3(ratio)))
  cat(sprintf("  J/I_etaeta = %.4g   Delta = %.3g\n",
              x$J_over_I_etaeta, x$Delta))
  if (is.finite(ratio) && abs(ratio - 1) > 0.1)
    cat("  Note: B differs from J/n by more than 10%; the expansion is local in gamma_V\n",
        "  and its leading term is not an accurate approximation at this value.\n", sep = "")
  invisible(x)
}

#' Plots for fitted bpss models
#'
#' \code{type = "score"} plots the efficient-score p-value against the
#' hypothesized correlation, with the nominal level as a dashed line and the
#' confidence set shaded; \code{type = "profile"} plots the relative profile
#' likelihood \eqn{R_n}{R} with reference levels 1/8 and 1/32;
#' \code{type = "information"} plots \eqn{\hat J_n / \hat I_{\eta\eta,n}}{J / I_etaeta}
#' and \eqn{\hat\Delta_n}{Delta} at the constrained fits. A dotted vertical line
#' marks the unrestricted estimate.
#'
#' @param x A fitted object of class \code{"bpss"}.
#' @param type Which plot to draw.
#' @param ... Further graphical parameters passed to \code{plot()}.
#' @return \code{x}, invisibly.
#' @export
plot.bpss <- function(x, type = c("score", "profile", "information"), ...) {
  type <- match.arg(type)
  gi <- x$grid
  if (is.null(gi)) stop("no grid inversion stored; refit with do_grid = TRUE.")
  rho_hat <- if (isTRUE(x$unrestricted$optim_ok) &&
                 isTRUE(x$unrestricted$likelihood_ok) &&
                 is.finite(x$unrestricted$eta_hat))
    tanh(x$unrestricted$eta_hat) else NA_real_
  xlab <- expression(rho^"*")
  if (type == "score") {
    curve <- bpss_score_curve(x)
    plot(curve$rho, curve$p, type = "l", xlab = xlab,
         ylab = "Score-test p-value", ylim = c(0, 1), ...)
    graphics::abline(h = 1 - x$level, lty = 2)
    set <- gi$score_rho_set
    if (!is.null(set) && nrow(set))
      for (i in seq_len(nrow(set)))
        graphics::rect(set$lower[i], 0, set$upper[i], 1 - x$level,
                       col = grDevices::adjustcolor("grey50", 0.35), border = NA)
    if (is.finite(rho_hat)) graphics::abline(v = rho_hat, lty = 3)
  } else if (type == "profile") {
    ps <- bpss_profile(x)
    plot(ps$curve$rho, ps$curve$R, type = "l", xlab = xlab,
         ylab = "Relative profile likelihood", ylim = c(0, 1), ...)
    graphics::abline(h = ps$levels, lty = 2)
    if (is.finite(rho_hat)) graphics::abline(v = rho_hat, lty = 3)
  } else {
    tab <- gi$table
    ok_ratio <- is.finite(tab$rho) & is.finite(tab$J) & is.finite(tab$I_etaeta) & tab$I_etaeta > 0
    ok_delta <- is.finite(tab$rho) & is.finite(tab$Delta)
    if (!any(ok_ratio) && !any(ok_delta))
      stop("the information diagnostics were not computed at any grid value.")
    oldpar <- graphics::par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
    on.exit(graphics::par(oldpar), add = TRUE)
    if (any(ok_ratio)) {
      plot(tab$rho[ok_ratio], tab$J[ok_ratio] / tab$I_etaeta[ok_ratio], type = "l",
           xlab = xlab, ylab = expression(hat(J)[n] / hat(I)[eta * eta * "," * n]), ...)
      if (is.finite(rho_hat)) graphics::abline(v = rho_hat, lty = 3)
    } else {
      graphics::plot.new(); graphics::title(main = "Information ratio not computed")
    }
    if (any(ok_delta)) {
      plot(tab$rho[ok_delta], tab$Delta[ok_delta], type = "l",
           xlab = xlab, ylab = expression(hat(Delta)[n]), ...)
      if (is.finite(rho_hat)) graphics::abline(v = rho_hat, lty = 3)
    } else {
      graphics::plot.new(); graphics::title(main = "Delta_n not computed")
    }
  }
  invisible(x)
}

#' Compare dependence inference across model specifications
#'
#' Tabulates, for several fitted specifications, the unrestricted estimate of
#' \eqn{\rho}{rho}, the Wald interval, the efficient-score confidence set and
#' the diagnostics of the test of \eqn{\rho = 0}{rho = 0}
#' (\eqn{\hat J_n(0)}{J(0)}, \eqn{\hat J_n(0) / \hat I_{\eta\eta,n}(0)}{J(0) / I_etaeta(0)},
#' \eqn{\hat\Delta_n(0)}{Delta(0)} and the score and likelihood-ratio
#' p-values), as in the application in the paper.
#'
#' @param fits A named list of fitted \code{"bpss"} objects.
#' @param digits Number of decimals for the reported sets.
#' @return A data frame with one row per specification.
#' @export
bpss_compare <- function(fits, digits = 3) {
  if (!is.list(fits) || is.null(names(fits)) || any(!nzchar(names(fits))) ||
      !all(vapply(fits, inherits, logical(1), "bpss")))
    stop("'fits' must be a named list of 'bpss' objects.")
  fmt <- function(set, computed = TRUE) {
    if (!isTRUE(computed)) return("not computed")
    if (is.null(set) || nrow(set) == 0L) return("empty")
    paste(sprintf("[%.*f, %.*f]", digits, set$lower,
                  digits, set$upper), collapse = " U ")
  }
  rows <- lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    ig <- f$ignorability
    audit_ok <- !is.null(f$unrestricted_audited) &&
                isTRUE(f$unrestricted_audited$audit_ok)
    data.frame(
      specification = nm,
      rho_hat = if (isTRUE(f$unrestricted$optim_ok) &&
                    isTRUE(f$unrestricted$likelihood_ok))
        round(f$unrestricted$rho_hat, digits) else NA_real_,
      multi_start = audit_ok,
      eta_near_bound = isTRUE(f$unrestricted$eta_boundary_hit),
      coef_near_bound = isTRUE(f$unrestricted$coef_boundary_hit),
      wald_rho = fmt(f$wald$rho, f$wald$usable),
      score_rho = if (!is.null(f$grid))
        fmt(f$grid$score_rho_set, TRUE) else "not computed",
      J0 = if (!is.null(ig) && isTRUE(ig$ok))
        signif(ig$J, 4) else NA_real_,
      J0_over_I0 = if (!is.null(ig) && isTRUE(ig$ok))
        signif(ig$J_over_I_etaeta, 3) else NA_real_,
      Delta0 = if (!is.null(ig) && isTRUE(ig$ok))
        signif(ig$Delta, 3) else NA_real_,
      p_score_0 = if (!is.null(ig) && isTRUE(ig$ok))
        signif(ig$p_chisq, 3) else NA_real_,
      p_lr_0 = if (!is.null(ig) && isTRUE(ig$lr_usable))
        signif(ig$p_lr, 3) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
