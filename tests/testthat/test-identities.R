# Algebraic identities and independent cross-checks at the reference design.

des <- paper_design(1L, 500L)
th <- des$theta

test_that("cell probabilities and their derivatives are coherent", {
  cd <- bpssInfo:::bpss_cell_derivatives_vec(th, des$XS, des$XO)
  expect_lt(max(abs(rowSums(cd$prob) - 1)), 1e-12)
  expect_lt(max(abs(apply(cd$D, c(1L, 3L), sum))), 1e-12)
})

test_that("compact and array information agree, and J is the Schur complement", {
  a <- bpssInfo:::bpss_info_array(th, des$XS, des$XO)
  b <- bpssInfo:::bpss_info_compact(th, des$XS, des$XO)
  expect_lt(max(abs(a$I - b$I)) / max(abs(a$I)), 1e-10)
  expect_lt(abs(a$J - b$J) / a$J, 1e-10)
  expect_lt(abs(b$J - b$J_schur) / b$J, 1e-8)
  expect_lt(max(abs(a$residualized_gradient - b$residualized_gradient)), 1e-10)
})

test_that("an independent block-form implementation gives the same information", {
  Inative <- bpss_expected_info_native(th, des$XS, des$XO)
  Ipkg <- bpssInfo:::bpss_info_compact(th, des$XS, des$XO)$I
  expect_lt(max(abs(Inative - unname(Ipkg))) / max(abs(Ipkg)), 1e-8)
})

test_that("the efficient score has mean zero and variance J under the model", {
  st <- bpssInfo:::bpss_efficient_score_stat(th, des$XS, des$XO,
                                          D = bpssInfo:::bpss_cell_prob(th, des$XS, des$XO))
  prob <- st$info$cell$prob
  expect_lt(abs(sum(prob * st$s)), 1e-10 * sqrt(st$J))
  expect_lt(abs(sum(prob * st$s^2) - st$J) / st$J, 1e-10)
})

test_that("the analytic exclusion-path derivative matches finite differences", {
  dec <- bpssInfo:::bpss_decomp_theta(th, des$XS, des$XO, idx_V = 3L,
                                  validate_fd = TRUE)
  expect_true(dec$fd_validation_pass)
  expect_equal(dec$B, dec$K0 + dec$cross_term + dec$exclusion_term)
  expect_equal(dec$A, dec$K0 + dec$exclusion_term)
})

test_that("a small selection probability is computed without cancellation", {
  pc <- bpssInfo:::.bpss_prob_core(a = c(0, 0), b = c(-6, -8), rho = 0.3)
  expect_lt(max(abs(pc$ps / stats::pnorm(c(-6, -8)) - 1)), 1e-12)
  expect_lt(max(abs(rowSums(pc$prob) - 1)), 1e-15)
})
