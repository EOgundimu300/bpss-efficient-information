# Values printed in Table 1 and Tables S4 and S5 of the paper (n = 2000 unless
# stated), reproduced at four significant figures from the designs of Section S5.1.

expect_signif <- function(x, printed, digits = 4) expect_equal(signif(x, digits), printed)

test_that("Tables 1 and S4 are reproduced at n = 2000", {
  S4 <- data.frame(
    J = c(13.14, 74.60, 0.8785, 4.801, 63.17, 5.597),
    ratio = c(0.08570, 0.4495, 0.01645, 0.03910, 0.5052, 0.02735),
    K0 = c(5.153e-3, 3.515e-5, 1.915e-5, 9.328e-10, 2.196e-6, 2.230e-3),
    cross = c(9.989e-5, 2.422e-4, 2.947e-7, 3.399e-9, 2.596e-5, 4.157e-6),
    excl = c(1.913e-3, 0.1636, 4.250e-4, 2.535e-3, 0.1341, 6.326e-4),
    kappa = c(0.01591, 0.05050, 0.001633, NA, 0.02392, 0.001750),
    omega = c(0.6093, 68.21, 4.711, NA, 247.1, 0.5327),
    B = c(7.165e-3, 0.1638, 4.445e-4, 2.535e-3, 0.1342, 2.866e-3),
    J_over_n = c(6.570e-3, 0.03730, 4.392e-4, 2.400e-3, 0.03159, 2.799e-3),
    d = c(-0.96, -0.22, 0.22, 0.02, -0.10, 0.82)
  )
  for (k in seq_len(nrow(paper_designs))) {
    des <- paper_design(k, 2000L)
    th <- des$theta
    beta <- th[4:5]; gamma <- th[1:3]; rho <- tanh(th[6])
    info <- bpss_information(beta, gamma, rho, des$XO, des$XS)
    expect_signif(info$J, S4$J[k])
    expect_signif(info$J_over_I_etaeta, S4$ratio[k])
    dec <- bpss_decomposition(beta, gamma, rho, des$XO, des$XS, exclusion = "V")
    expect_signif(dec$K0, S4$K0[k])
    expect_signif(dec$cross_term, S4$cross[k])
    expect_signif(dec$exclusion_term, S4$excl[k])
    expect_signif(dec$B, S4$B[k])
    expect_signif(dec$J_over_n, S4$J_over_n[k])
    if (is.na(S4$kappa[k])) {
      expect_true(is.na(dec$kappa) && is.na(dec$omega))
    } else {
      expect_signif(dec$kappa, S4$kappa[k])
      expect_signif(dec$omega, S4$omega[k])
    }
    p <- paper_designs[k, ]
    expect_equal(round(p$gammax + des$corr_emp * p$gammaV - p$rho0 * p$beta1, 2), S4$d[k])
  }
})

test_that("Table S5 is reproduced at n = 500, 1000 and 4000", {
  S5 <- rbind(c(3.294, 6.576, 26.27), c(18.71, 37.34, 149.1), c(0.2206, 0.4399, 1.756),
              c(1.204, 2.403, 9.595), c(15.84, 31.63, 126.3), c(1.400, 2.799, 11.19))
  for (k in seq_len(nrow(paper_designs))) for (j in 1:3) {
    des <- paper_design(k, c(500L, 1000L, 4000L)[j])
    th <- des$theta
    info <- bpss_information(th[4:5], th[1:3], tanh(th[6]), des$XO, des$XS)
    expect_signif(info$J, S5[k, j])
  }
})
