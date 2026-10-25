# Independent block-form expected information, used only to check the
# package's implementation. Parameter order: (gamma, beta, eta).
bpss_expected_info_native <- function(beta, XS, XO) {
  XS <- as.matrix(XS); XO <- as.matrix(XO)
  NXS <- ncol(XS); NXO <- ncol(XO); n <- nrow(XS)
  ibetaS <- seq_len(NXS); ibetaO <- NXS + seq_len(NXO); irho <- NXS + NXO + 1L; nParam <- irho
  gammaS <- beta[ibetaS]; betaO <- beta[ibetaO]; rho <- tanh(beta[irho]); r <- sqrt(1 - rho^2)
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  mu <- drop(XO %*% betaO); eta <- drop(XS %*% gammaS)
  omg1 <- (eta - rho * mu) / r; omg2 <- (-mu + rho * eta) / r
  phi_mu <- stats::dnorm(mu); phi_eta <- stats::dnorm(eta)
  Phi_omg1 <- stats::pnorm(omg1); Phi_omg2 <- stats::pnorm(omg2); Phi_nomg2 <- stats::pnorm(-omg2)
  p0 <- pmax(stats::pnorm(-eta), 1e-16); p11 <- pmax(as.numeric(pbivnorm::pbivnorm(eta, mu, rho)), 1e-16)
  p10 <- pmax(stats::pnorm(eta) - p11, 1e-16); phi2 <- mvtnorm::dmvnorm(cbind(eta, mu), sigma = Sigma)
  A <- phi_mu * Phi_omg1; B <- phi_eta * Phi_nomg2; C <- r^2 * phi2; Dd <- phi_eta * Phi_omg2
  w_bb <- A^2 * (1 / p10 + 1 / p11); w_gg <- phi_eta^2 / p0 + Dd^2 / p10 + B^2 / p11; w_aa <- C^2 * (1 / p10 + 1 / p11)
  w_bg <- phi_mu * phi_eta * Phi_omg1 * (Phi_nomg2 / p11 - Phi_omg2 / p10)
  w_ba <- A * C * (1 / p10 + 1 / p11); w_ga <- phi_eta * C * (Phi_nomg2 / p11 - Phi_omg2 / p10)
  FI <- matrix(0, nParam, nParam)
  FI[ibetaO, ibetaO] <- crossprod(XO * as.vector(w_bb), XO)
  FI[ibetaS, ibetaS] <- crossprod(XS * as.vector(w_gg), XS); FI[irho, irho] <- sum(w_aa)
  FI[ibetaO, ibetaS] <- crossprod(XO * as.vector(w_bg), XS); FI[ibetaS, ibetaO] <- t(FI[ibetaO, ibetaS])
  FI[ibetaO, irho] <- drop(crossprod(XO, w_ba)); FI[irho, ibetaO] <- FI[ibetaO, irho]
  FI[ibetaS, irho] <- drop(crossprod(XS, w_ga)); FI[irho, ibetaS] <- FI[ibetaS, irho]
  FI
}
