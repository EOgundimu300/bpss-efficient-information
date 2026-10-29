# User-facing fitting function and numerical controls.

# Coerce a binary response to 0/1 integer, preserving NA.
.bpss_binary <- function(z, what = "response") {
  if (is.logical(z)) return(as.integer(z))
  if (is.factor(z)) {
    if (nlevels(z) != 2L) stop(what, " must have exactly two levels.")
    return(as.integer(z) - 1L)
  }
  zn <- suppressWarnings(as.numeric(z))
  u <- sort(unique(zn[is.finite(zn)]))
  if (length(u) && all(u %in% c(0, 1))) return(as.integer(zn))
  stop(what, " must be 0/1, logical, or a two-level factor.")
}

#' Bivariate probit sample selection model with efficient-score inference
#'
#' Fits the bivariate probit sample selection model and inverts the efficient
#' score test of \eqn{\eta = \mathrm{atanh}(\rho)}{eta = atanh(rho)} on a grid of
#' hypothesized correlations.
#'
#' The model has latent equations
#' \eqn{Y^* = x^\top\beta + \varepsilon_1}{Y* = x'beta + e1} for the outcome and
#' \eqn{S^* = w^\top\gamma + \varepsilon_2}{S* = w'gamma + e2} for selection, where
#' \eqn{(\varepsilon_1, \varepsilon_2)}{(e1, e2)} is standard bivariate normal
#' with correlation \eqn{\rho}{rho}. The selection indicator
#' \eqn{S = 1(S^* > 0)}{S = 1(S* > 0)} is observed for every unit, and the outcome
#' \eqn{Y = 1(Y^* > 0)}{Y = 1(Y* > 0)} only when \eqn{S = 1}{S = 1}.
#'
#' At each hypothesized value \eqn{\eta^\star}{eta*} the regression coefficients
#' are estimated with \eqn{\eta}{eta} fixed, and
#' \eqn{T_n(\eta^\star) = \hat U_n(\eta^\star)^2 / \hat J_n(\eta^\star)}{T(eta*) = U(eta*)^2 / J(eta*)}
#' is referred to the \eqn{\chi^2_1}{chi-squared(1)} distribution, where
#' \eqn{\hat U_n}{U} is the efficient score and \eqn{\hat J_n}{J} the expected
#' efficient information at the constrained fit. The confidence set consists of
#' the evaluated grid values at which the test does not reject, with endpoints
#' refined where the decision changes between neighbouring grid values. A
#' component narrower than the grid spacing can be missed; a finer grid
#' (\code{grid_length} in \code{\link{bpss_control}}, or \code{eta_grid})
#' reduces this risk. Grid values at which the statistic cannot be computed are
#' reported, and a component reaching an end of the grid is reported only up to
#' that end. The paper recommends reporting, for each component that does not
#' contain the estimate or that reaches an end of the grid, the relative profile
#' likelihood and the ranges of \eqn{\hat J_n}{J},
#' \eqn{\hat J_n / \hat I_{\eta\eta,n}}{J / I_etaeta} and
#' \eqn{\hat\Delta_n}{Delta}; \code{summary()} prints these.
#'
#' The unrestricted maximum-likelihood fit is obtained by bounded L-BFGS-B
#' optimization and, by default, a multi-start search over fixed-correlation
#' starts. It is the centre of the grid continuation, the reference for the
#' likelihood-ratio statistics and the relative profile likelihood, and the
#' centre of the Wald interval.
#'
#' The efficient-score set is the primary inference on \eqn{\rho}{rho}. The
#' Wald interval relies on a quadratic approximation to the log-likelihood at
#' the unrestricted estimate, can be unreliable when the efficient information
#' is small, and is reported for comparison. It is computed only when the
#' unrestricted fit passes its convergence, stationarity and information checks,
#' no parameter is near its bound and no constrained fit has a larger
#' log-likelihood; the further profile-likelihood condition applied to Wald
#' intervals in the simulations of the paper is not implemented.
#'
#' Covariates must be observed for every unit, including unselected units,
#' because the expected information sums over all units.
#'
#' @param selection Formula for the selection equation, with the selection
#'   indicator (0/1, logical or a two-level factor) on the left-hand side.
#' @param outcome Formula for the outcome equation, with the binary outcome
#'   (0/1, logical or a two-level factor) on the left-hand side; the outcome may
#'   be \code{NA} when the selection indicator is 0.
#' @param data A data frame.
#' @param subset Optional expression selecting rows of \code{data}.
#' @param eta_grid Optional grid of \eqn{\eta}{eta} values; by default
#'   \code{control$grid_length} equally spaced values of \eqn{\rho}{rho} on
#'   \code{control$grid_rho}. The value 0 is always added.
#' @param level Confidence level.
#' @param control Numerical controls from \code{\link{bpss_control}}.
#'
#' @return An object of class \code{"bpss"}, a list including the unrestricted
#'   fits (\code{unrestricted}, \code{unrestricted_default},
#'   \code{unrestricted_audited}), the grid inversion (\code{grid}, whose
#'   \code{table} holds at each grid value the statistic, its p-value,
#'   \eqn{\hat J_n}{J}, \eqn{\hat I_{\eta\eta,n}}{I_etaeta}, \eqn{\hat\Delta_n}{Delta},
#'   the smallest fitted cell probability and \eqn{-2\log R_n}{-2 log R}), the
#'   Wald interval (\code{wald}), the test of \eqn{\rho = 0}{rho = 0}
#'   (\code{ignorability}) and the design matrices.
#'
#' @references Ogundimu, E. O. Efficient information for the dependence
#'   parameter in bivariate probit sample selection models. Manuscript.
#' @seealso \code{\link{bpss_control}}, \code{\link{summary.bpss}},
#'   \code{\link{confint.bpss}}, \code{\link{bpss_score_test}},
#'   \code{\link{test_ignorability}}, \code{\link{bpss_profile}},
#'   \code{\link{bpss_decompose}}
#' @examples
#' set.seed(1)
#' n <- 800
#' x <- rnorm(n)                            # covariate in both equations
#' v <- rnorm(n)                            # selection-only covariate
#' e1 <- rnorm(n)
#' e2 <- 0.5 * e1 + sqrt(0.75) * rnorm(n)   # errors with correlation 0.5
#' s <- as.integer(0.3 + 0.5 * x + 0.8 * v + e2 > 0)
#' y <- ifelse(s == 1, as.integer(-0.4 + 0.6 * x + e1 > 0), NA)
#' dat <- data.frame(s, y, x, v)
#' fit <- bpss(s ~ x + v, y ~ x, data = dat)
#' summary(fit)
#' confint(fit)
#' plot(fit, "score")
#' @export
bpss <- function(selection, outcome, data, subset = NULL,
                 eta_grid = NULL, level = 0.95,
                 control = bpss_control()) {
  cl <- match.call()
  if (!inherits(selection, "formula") ||
      !inherits(outcome, "formula"))
    stop("'selection' and 'outcome' must be formulas.")
  if (missing(data) || !is.data.frame(data))
    stop("'data' must be supplied as a data.frame.")

  subset_expr <- match.call(expand.dots = FALSE)$subset
  if (!is.null(subset_expr)) {
    row_idx <- eval(subset_expr, data, parent.frame())
    data <- data[row_idx, , drop = FALSE]
  }
  mfS <- stats::model.frame(
    selection, data = data, na.action = stats::na.pass,
    drop.unused.levels = TRUE
  )
  mfO <- stats::model.frame(
    outcome, data = data, na.action = stats::na.pass,
    drop.unused.levels = TRUE
  )
  if (nrow(mfS) != nrow(mfO))
    stop("selection and outcome frames have different row counts.")

  S_raw <- stats::model.response(mfS)
  Y_raw <- stats::model.response(mfO)
  if (is.null(S_raw))
    stop("selection formula must have a left-hand side.")
  if (is.null(Y_raw))
    stop("outcome formula must have a left-hand side.")
  S01 <- .bpss_binary(S_raw, "selection indicator")

  predS <- if (ncol(mfS) > 1L)
    mfS[, -1L, drop = FALSE] else mfS[, 0L, drop = FALSE]
  predO <- if (ncol(mfO) > 1L)
    mfO[, -1L, drop = FALSE] else mfO[, 0L, drop = FALSE]
  complete_pred <- stats::complete.cases(predS) &
                   stats::complete.cases(predO)
  selected <- !is.na(S01) & S01 == 1L
  if (any(is.na(S01)))
    warning(sum(is.na(S01)),
            " row(s) dropped for a missing selection indicator.")
  drop_cov <- !is.na(S01) & !complete_pred
  if (any(drop_cov))
    warning(sum(drop_cov), " row(s) dropped for missing covariates (",
            sum(S01[drop_cov] == 1L), " selected, ",
            sum(S01[drop_cov] == 0L), " unselected); outcome covariates ",
            "are also needed for unselected rows.")
  drop_selna <- selected & is.na(Y_raw) & complete_pred
  if (any(drop_selna))
    warning(sum(drop_selna),
            " selected row(s) dropped for missing outcome.")
  keep <- complete_pred & !is.na(S01) &
          !(selected & is.na(Y_raw))
  if (sum(keep) < 10L)
    stop("fewer than 10 usable rows after missing-data handling.")

  XS <- stats::model.matrix(
    attr(mfS, "terms"), mfS[keep, , drop = FALSE])
  XO <- stats::model.matrix(
    attr(mfO, "terms"), mfO[keep, , drop = FALSE])
  S01 <- S01[keep]
  Y01 <- .bpss_binary(Y_raw[keep], "outcome")
  Y01[S01 == 0L] <- 0L
  D <- observed_cell_matrix(S01, Y01)
  dd <- theta_dims(XS, XO)

  nm_S <- setdiff(colnames(XS), "(Intercept)")
  nm_O <- setdiff(colnames(XO), "(Intercept)")
  excluded <- setdiff(nm_S, nm_O)

  fit_default <- bpss_fit_unrestricted(
    D, XS, XO, eps = control$eps,
    rcond_tol = control$rcond_tol,
    coef_bound = control$coef_bound,
    eta_bound = control$eta_bound,
    maxit = control$maxit,
    optim_factr = control$optim_factr,
    optim_pgtol = control$optim_pgtol,
    stationarity_tol = control$stationarity_tol,
    boundary_fraction = control$boundary_fraction
  )
  fit_audited <- NULL
  if (isTRUE(control$audit_unrestricted)) {
    fit_audited <- bpss_fit_unrestricted_audited(
      D, XS, XO, default_fit = fit_default,
      eta_starts = atanh(control$audit_eta_rho),
      use_constrained_starts =
        control$audit_use_constrained_starts,
      eps = control$eps, rcond_tol = control$rcond_tol,
      coef_bound = control$coef_bound,
      eta_bound = control$eta_bound,
      maxit = control$maxit,
      optim_factr = control$optim_factr,
      optim_pgtol = control$optim_pgtol,
      stationarity_tol = control$stationarity_tol,
      boundary_fraction = control$boundary_fraction,
      ll_tol = control$audit_ll_tol,
      par_tol = control$audit_par_tol
    )
  }
  fitU <- if (!is.null(fit_audited)) fit_audited else fit_default

  if (is.null(eta_grid)) {
    rho_grid <- seq(
      control$grid_rho[1], control$grid_rho[2],
      length.out = control$grid_length
    )
    eta_grid <- atanh(rho_grid)
  }
  eta_grid <- sort(unique(c(as.numeric(eta_grid), 0)))

  # The profile likelihood uses fitU as reference when it passed the audit.
  grid <- NULL
  if (isTRUE(control$do_grid)) {
    grid <- bpss_grid_invert(
      D, XS, XO, eta_grid, unrestricted = fitU,
      level = level,
      eps = control$eps, rcond_tol = control$rcond_tol,
      coef_bound = control$coef_bound,
      maxit = control$maxit,
      optim_factr = control$optim_factr,
      optim_pgtol = control$optim_pgtol,
      stationarity_tol = control$stationarity_tol,
      boundary_fraction = control$boundary_fraction,
      refine = control$refine_grid,
      refine_tol = control$refine_tol,
      keep_fits = TRUE,
      constrained_audit = control$audit_constrained_grid,
      audit_ll_tol = control$audit_ll_tol,
      audit_par_tol = control$audit_par_tol,
      require_audited_profile =
        control$require_audited_profile,
      profile_lr_tol = control$profile_lr_tol
    )
  }

  # A constrained fit with a larger log-likelihood shows that fitU is not the
  # maximum; it then supports neither Wald nor likelihood-ratio inference.
  if (isTRUE(grid$profile_reference_dominated)) {
    fitU$audit_ok <- FALSE
    fitU$reference_dominated <- TRUE
  }

  obj <- list(
    call = cl, selection = selection, outcome = outcome,
    terms = list(
      selection = attr(mfS, "terms"),
      outcome = attr(mfO, "terms")
    ),
    xlevels = list(
      selection = stats::.getXlevels(attr(mfS, "terms"), mfS),
      outcome = stats::.getXlevels(attr(mfO, "terms"), mfO)
    ),
    XS = XS, XO = XO, D = D, S = S01, Y = Y01,
    dims = dd, coef_names = make_theta_names(XS, XO),
    excluded = excluded,
    control = control, level = level,
    unrestricted_default = fit_default,
    unrestricted_audited = fit_audited,
    unrestricted = fitU, grid = grid
  )
  class(obj) <- "bpss"
  obj$wald <- bpss_wald(obj, level = level)
  obj$ignorability <- tryCatch(
    test_ignorability(obj, require_audited_lr = TRUE),
    error = function(e)
      list(ok = FALSE, status = conditionMessage(e))
  )
  obj
}

#' Numerical controls for bpss
#'
#' @param eps Smallest admissible fitted cell probability for score and Wald
#'   statistics.
#' @param rcond_tol Smallest admissible reciprocal condition number of the
#'   nuisance information.
#' @param coef_bound Bound on the absolute value of each regression
#'   coefficient during optimization.
#' @param eta_bound Bound on \eqn{|\eta|}{|eta|}; the default corresponds to
#'   \eqn{|\rho| \le 0.995}{|rho| <= 0.995}.
#' @param maxit,optim_factr,optim_pgtol Settings of the first L-BFGS-B run; a
#'   failed or nonstationary run is retried with \code{2 * maxit} iterations
#'   and \code{factr = min(optim_factr, 1e4)}.
#' @param stationarity_tol Largest admissible absolute projected-gradient
#'   component of the negative average log-likelihood.
#' @param boundary_fraction A parameter is classed as near its bound when its
#'   absolute value is at least this fraction of the bound.
#' @param grid_length,grid_rho Number and range of equally spaced values of
#'   \eqn{\rho}{rho} in the inversion grid.
#' @param refine_grid,refine_tol Refine set endpoints at changes in the
#'   decision, with this root-finding tolerance on the \eqn{\eta}{eta} scale.
#' @param do_grid Logical; invert the score test on the grid.
#' @param audit_unrestricted Logical; obtain the unrestricted fit by
#'   multi-start search. Required for the likelihood-ratio statistics and the
#'   relative profile likelihood.
#' @param audit_eta_rho Correlations at which constrained fits supply starting
#'   values for the multi-start search.
#' @param audit_use_constrained_starts Logical; use those constrained fits as
#'   starts (otherwise separate-probit starts at the same correlations).
#' @param audit_ll_tol,audit_ll_rel_tol Two log-likelihoods are treated as
#'   different when they differ by more than
#'   \code{max(audit_ll_tol, audit_ll_rel_tol * max(1, |logLik|))}.
#' @param audit_rho_tol,audit_eta_se_tol Thresholds on differences in
#'   \eqn{\hat\rho}{rho-hat}, and in \eqn{\hat\eta}{eta-hat} in standard-error units, that
#'   trigger the remaining starts and the discordance flags.
#' @param audit_top_starts Number of constrained fits with the largest
#'   log-likelihoods used as initial unrestricted starts, in addition to the
#'   single-start fit and the two extreme correlations.
#' @param audit_expand Logical; add the remaining starts when the initial
#'   unrestricted solutions disagree or are exceeded by a constrained fit.
#' @param audit_par_tol Parameter difference above which fits with equal
#'   log-likelihoods are flagged as distinct.
#' @param audit_constrained_grid Audit of constrained fits on the grid:
#'   \code{"adaptive"} (default), \code{"all"} or \code{"none"}.
#' @param grid_audit_stride Spacing of the regular audit checkpoints.
#' @param require_audited_profile Logical; compute likelihood-ratio
#'   statistics only against an audited unrestricted fit.
#' @param profile_lr_tol A constrained fit whose log-likelihood exceeds that of
#'   the reference fit by more than half this value invalidates the profile
#'   likelihood.
#'
#' @return A list of control values.
#' @export
bpss_control <- function(
    eps = 1e-12, rcond_tol = 1e-12, coef_bound = 12,
    eta_bound = atanh(0.995), maxit = 1000,
    optim_factr = 1e7, optim_pgtol = 1e-8,
    stationarity_tol = 1e-4, boundary_fraction = 0.99,
    grid_length = 201L, grid_rho = c(-0.99, 0.99),
    refine_grid = TRUE, refine_tol = 1e-6, do_grid = TRUE,
    audit_unrestricted = TRUE,
    audit_eta_rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    audit_use_constrained_starts = TRUE,
    audit_ll_tol = 1e-3, audit_ll_rel_tol = 1e-8,
    audit_rho_tol = 0.01, audit_eta_se_tol = 0.05,
    audit_top_starts = 3L, audit_expand = TRUE,
    audit_par_tol = 1e-3,
    audit_constrained_grid = "adaptive", grid_audit_stride = 20L,
    require_audited_profile = TRUE, profile_lr_tol = 1e-3) {
  if (is.logical(audit_constrained_grid)) {
    audit_mode <- if (isTRUE(audit_constrained_grid)) "all" else "none"
  } else {
    audit_mode <- match.arg(as.character(audit_constrained_grid)[1L],
                            c("adaptive", "all", "none"))
  }
  stopifnot(
    length(grid_rho) == 2L, grid_rho[1] < grid_rho[2],
    all(abs(grid_rho) < 1), grid_length >= 11L,
    all(abs(audit_eta_rho) < 1), refine_tol > 0,
    optim_pgtol >= 0, stationarity_tol > 0,
    boundary_fraction > 0, boundary_fraction < 1,
    audit_ll_tol > 0, audit_ll_rel_tol >= 0,
    audit_rho_tol > 0, audit_eta_se_tol > 0,
    audit_top_starts >= 1L, audit_par_tol > 0,
    grid_audit_stride >= 1L, profile_lr_tol >= 0
  )
  ll_control <- structure(
    as.numeric(audit_ll_tol),
    rel_tol = audit_ll_rel_tol,
    rho_tol = audit_rho_tol,
    eta_se_tol = audit_eta_se_tol,
    top_starts = as.integer(audit_top_starts),
    expand = isTRUE(audit_expand)
  )
  grid_audit <- structure(audit_mode,
                          stride = as.integer(grid_audit_stride))
  list(
    eps = eps, rcond_tol = rcond_tol,
    coef_bound = coef_bound, eta_bound = eta_bound,
    maxit = maxit, optim_factr = optim_factr,
    optim_pgtol = optim_pgtol,
    stationarity_tol = stationarity_tol,
    boundary_fraction = boundary_fraction,
    grid_length = as.integer(grid_length), grid_rho = grid_rho,
    refine_grid = refine_grid, refine_tol = refine_tol,
    do_grid = do_grid, audit_unrestricted = audit_unrestricted,
    audit_eta_rho = audit_eta_rho,
    audit_use_constrained_starts = audit_use_constrained_starts,
    audit_ll_tol = ll_control, audit_par_tol = audit_par_tol,
    audit_constrained_grid = grid_audit,
    require_audited_profile = require_audited_profile,
    profile_lr_tol = profile_lr_tol
  )
}
