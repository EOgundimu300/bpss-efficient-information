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

## Example

The selection and outcome equations are given as formulas. The outcome may be
`NA` when the selection indicator is 0; covariates must be observed for all
units. The example simulates data with error correlation 0.5 and a covariate
`v` that enters only the selection equation.

```r
library(bpssInfo)

set.seed(1)
n  <- 800
x  <- rnorm(n)                          # covariate in both equations
v  <- rnorm(n)                          # selection-only covariate
e1 <- rnorm(n)
e2 <- 0.5 * e1 + sqrt(0.75) * rnorm(n)  # errors with correlation 0.5
s  <- as.integer(0.3 + 0.5 * x + 0.8 * v + e2 > 0)              # selection indicator
y  <- ifelse(s == 1, as.integer(-0.4 + 0.6 * x + e1 > 0), NA)    # outcome, observed if s = 1
dat <- data.frame(s, y, x, v)

fit <- bpss(s ~ x + v, y ~ x, data = dat)
summary(fit)          # estimates, information diagnostics, score set and its components
confint(fit)          # efficient-score confidence set for rho (default)
plot(fit, "score")    # score-test p-value against the hypothesized correlation
plot(fit, "profile")  # relative profile likelihood
test_ignorability(fit)            # score and likelihood-ratio tests of rho = 0
bpss_score_test(fit, rho = 0.3)   # score test at a specified correlation
```

`bpss_information()` and `bpss_decomposition()` evaluate the efficient
information and its exclusion decomposition at specified coefficients and
design matrices, as in the design-level calculations of the paper;
`bpss_decompose()` evaluates the decomposition for a fitted model. The
decomposition is a local expansion in the coefficient of the excluded
covariate, and its leading term need not approximate the information when that
coefficient is not small; the printed ratio B/(J/n) shows how close it is.
`bpss_compare()` tabulates the results for several model specifications.

The score confidence set is computed on a grid of correlations (201 values on
[-0.99, 0.99] by default), with endpoints refined where the decision changes.
Grid values at which the statistic cannot be computed are reported, and a
finer grid can be requested with `bpss_control(grid_length = )`. The
efficient-score set is the primary inference on the correlation. The
unrestricted fit uses a multi-start search, and the Wald interval centred at it
is reported for comparison; it can be unreliable when the efficient information
is small. The Wald interval is computed only when that fit passes its checks,
no parameter is near its bound and no constrained fit has a larger
log-likelihood; the further profile-likelihood condition applied to Wald
intervals in the simulations of the paper is not implemented.

The tests in `tests/testthat` reproduce design-level information values
reported in the paper and check the information identities, the zero-information
locus and fourth-order decay of the efficient information, the analytic
derivative used in the exclusion decomposition and the factorized fit at
independence.

## Licence

GPL (>= 2).
