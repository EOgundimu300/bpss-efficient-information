# Local exclusion decomposition and design-level information.

test_that("the leading term approximates J/n for a small exclusion coefficient", {
  n <- 500
  x <- seq(-1.5, 1.5, length.out = n)
  XO <- cbind("(Intercept)" = 1, x = x)
  XS <- cbind("(Intercept)" = 1, x = x, v = x^2 - mean(x^2))
  dec <- bpss_decomposition(c(-0.5, 0.4), c(0.4, 0.6, 0.05), 0.3, XO, XS,
                            exclusion = "v")
  expect_lt(abs(dec$B / dec$J_over_n - 1), 0.01)
  expect_lt(dec$kappa, 0)
  expect_output(print(dec), "negative local cross term")
})

test_that("bpss_information uses the exact probabilities when eps = 0", {
  des <- paper_design(1L, 2000L)
  g0 <- .safe_uniroot(function(g) mean(stats::pnorm(g + 1.5 * des$x)) - 0.65)
  beta <- c(des$params$beta0, des$params$beta1)
  gamma <- c(g0, 1.5, 0)
  expect_error(bpss_information(beta, gamma, -0.9, des$XO, des$XS), "eps = 0")
  exact <- bpss_information(beta, gamma, -0.9, des$XO, des$XS, eps = 0)
  floored <- bpssInfo:::bpss_info_compact(c(gamma, beta, atanh(-0.9)), des$XS, des$XO)
  expect_lt(abs(exact$J_over_I_etaeta - floored$J / floored$I_etaeta), 1e-8)
})
