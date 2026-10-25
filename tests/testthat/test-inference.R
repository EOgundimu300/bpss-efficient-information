# Consistency of the inferential output on simulated data.

set.seed(20260713)
n <- 600
x <- rnorm(n); v <- rnorm(n)
e1 <- rnorm(n); e2 <- 0.5 * e1 + sqrt(0.75) * rnorm(n)
s <- as.integer(0.3 + 0.5 * x + 0.8 * v + e2 > 0)
y <- ifelse(s == 1, as.integer(-0.4 + 0.6 * x + e1 > 0), NA)
dat <- data.frame(s, y, x, v)
fit <- bpss(s ~ x + v, y ~ x, data = dat,
            control = bpss_control(grid_length = 31))

test_that("bpss returns a fitted object with a computed score set", {
  expect_s3_class(fit, "bpss")
  expect_true(isTRUE(fit$unrestricted$audit_ok))
  ci <- confint(fit)
  expect_s3_class(ci, "bpss_set")
  expect_gt(nrow(ci), 0L)
  expect_true(all(ci$lower < ci$upper))
  expect_output(print(summary(fit)), "efficient-score set")
})

test_that("the score test at a grid value agrees with the grid inversion", {
  tab <- fit$grid$table
  j <- which(tab$resolved & abs(tab$rho - 0.4) < 0.05)[1L]
  st <- bpss_score_test(fit, tab$rho[j])
  expect_true(st$computed)
  expect_equal(st$T, tab$T[j], tolerance = 1e-4)
  expect_equal(st$J, tab$J[j], tolerance = 1e-4)
})

test_that("the test of rho = 0 agrees with the score test at zero", {
  ig <- test_ignorability(fit)
  st <- bpss_score_test(fit, 0)
  expect_true(ig$ok)
  expect_equal(st$T, ig$T, tolerance = 1e-4)
  expect_equal(st$p_value, ig$p_chisq, tolerance = 1e-4)
})

test_that("the relative profile likelihood is at most one", {
  pr <- bpss_profile(fit)
  R <- pr$curve$R[is.finite(pr$curve$R)]
  expect_true(all(R <= 1 + 1e-12))
  expect_gt(max(R), 0.9)
})

test_that("missing covariates are reported", {
  dat2 <- dat; dat2$x[c(1, 2)] <- NA
  expect_warning(bpss(s ~ x + v, y ~ x, data = dat2,
                      control = bpss_control(do_grid = FALSE,
                                             audit_unrestricted = FALSE)),
                 "missing covariates")
})

test_that("the factorized fit at rho = 0 is the constrained maximum", {
  nf <- bpssInfo:::bpss_fit_null_factorized(fit$D, fit$XS, fit$XO,
                                             check_generic = TRUE)
  expect_true(nf$ok)
  expect_identical(nf$status, "ok")
  expect_true(isTRUE(nf$generic_agrees))
  expect_lt(abs(nf$factorization_error), 1e-8)
  gc <- bpssInfo:::bpss_fit_constrained(0, fit$D, fit$XS, fit$XO)
  expect_lt(abs(gc$logLik - nf$logLik), 1e-6)
  expect_equal(gc$T_score, nf$T_score, tolerance = 1e-4)
})

test_that("the smallest cell probability is recorded when the statistic is not computed", {
  fc <- bpssInfo:::bpss_fit_constrained(atanh(0.99), fit$D, fit$XS, fit$XO)
  expect_true(is.finite(fc$min_cell_probability))
  if (isTRUE(grepl("probability", fc$information_error)))
    expect_lt(fc$min_cell_probability, 1e-12)
})

test_that("the likelihood-ratio test of rho = 0 is computed when eta is near its bound", {
  f2 <- fit
  f2$unrestricted$eta_boundary_hit <- TRUE
  ig2 <- test_ignorability(f2)
  expect_true(ig2$lr_usable)
  expect_true(ig2$lr_reference_at_bound)
  expect_equal(ig2$LR, test_ignorability(fit)$LR)
  f2$ignorability <- ig2
  expect_output(print(summary(f2)), "near the bound")
})
