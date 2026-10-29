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

test_that("a cell probability that is zero in double precision gives an informative error", {
  x <- seq(-3, 3, length.out = 200)
  XO <- cbind("(Intercept)" = 1, x = x)
  XS <- cbind("(Intercept)" = 1, x = x, v = sin(2 * x))
  expect_error(bpss_information(c(0, 1), c(0, 15, 0.3), 0.5, XO, XS, eps = 0),
               "zero in double precision")
})

test_that("the efficient information is zero on the locus gamma_x = rho beta_1", {
  des <- paper_design(1L, 500L)
  th <- des$theta
  rho0 <- tanh(th[6])
  info <- bpss_information(th[4:5], c(th[1], rho0 * th[5], 0), rho0, des$XO, des$XS)
  expect_lt(info$J_over_I_etaeta, 1e-12)
})

test_that("the efficient information decays as the fourth power of the slope departure", {
  des <- paper_design(1L, 500L)
  th <- des$theta
  rho0 <- tanh(th[6])
  J <- function(d)
    bpss_information(th[4:5], c(th[1], rho0 * th[5] + d, 0), rho0, des$XO, des$XS)$J
  expect_equal(J(0.02) / J(0.01), 16, tolerance = 0.01)
})

test_that("the decomposition requires cell probabilities above eps at the evaluation point", {
  x <- seq(-1, 1, length.out = 200)
  XO <- cbind("(Intercept)" = 1, x = x)
  XS <- cbind("(Intercept)" = 1, x = x, v = 3 * sin(3 * x))
  expect_error(bpss_decomposition(c(0, 0.5), c(0, 0.5, 4), 0.3, XO, XS,
                                  exclusion = "v"),
               "evaluation point")
})
