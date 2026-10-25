# bpss-efficient-information

R package `bpssInfo`, accompanying

> Ogundimu, E. O. Efficient information for the dependence parameter in
> bivariate probit sample selection models.

The package fits the bivariate probit sample selection model and makes
inference on the error correlation by inverting the efficient score test, with
the regression coefficients estimated at each hypothesized correlation. It
computes the efficient information for the Fisher-transformed correlation, its
local decomposition under an exclusion restriction, the relative profile
likelihood and the single-observation diagnostics recommended in the paper.

## Installation

```r
# install.packages("remotes")
remotes::install_github("EOgundimu300/bpss-efficient-information")
```

Requirements: R (>= 4.0.0) and the packages `pbivnorm` and `mvtnorm`.

## Use

The selection and outcome equations are given as formulas; the outcome may be
`NA` when the selection indicator is 0. Covariates must be observed for all
units.

```r
library(bpssInfo)

fit <- bpss(s ~ x + v, y ~ x, data = dat)
summary(fit)          # estimates, information diagnostics, score set and its components
confint(fit)          # efficient-score confidence set for rho (default)
plot(fit, "score")    # score-test p-value against the hypothesized correlation
plot(fit, "profile")  # relative profile likelihood
test_ignorability(fit)            # score and likelihood-ratio tests of rho = 0
bpss_score_test(fit, rho = 0.3)   # score test at a specified correlation
bpss_decompose(fit)               # local exclusion decomposition
```

`bpss_information()` and `bpss_decomposition()` evaluate the efficient
information and its exclusion decomposition at specified coefficients and
design matrices, as in the design-level calculations of the paper.
`bpss_compare()` tabulates the results for several model specifications.

The unrestricted fit uses the multi-start search of Section S5.1 of the paper.
The Wald interval is computed only when that fit passes its checks and no
parameter is near its bound; the further profile-likelihood condition applied
to Wald intervals in the simulations of the paper is not implemented.

The tests in `tests/testthat` reproduce Table 1 and Tables S4 and S5 of the
paper and check the information identities and the analytic derivative used
in the exclusion decomposition.

## Licence

GPL (>= 2).
