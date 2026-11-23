# QEAkash and QR-QEAkash with GAMLSS

`QEAkash_GAMLSS.R` defines the Quantile Exponentiated Akash (QEAkash)
family for the R package `gamlss`. Use an intercept-only model to fit a
common quantile and `sigma`, or include covariates in either parameter to
fit a QR-QEAkash regression.

## Setup

Install the required packages once:

```r
install.packages(c("gamlss", "gamlss.dist"))
```

Run the examples from the directory containing both R scripts:

```sh
Rscript QEAkash_Examples.R
```

To work interactively from that directory:

```r
library(gamlss)
source("QEAkash_GAMLSS.R")
```

## Family and parameterization

The response must be finite and positive. The Akash distribution with
parameter `sigma > 0` has distribution function

$$
F(x;\sigma) = 1 - \left(1 + \frac{\sigma x(\sigma x + 2)}{\sigma^2 + 2}\right)e^{-\sigma x},
\qquad x > 0.
$$

The QEAkash distribution raises it to the power
$\alpha = \log(\tau)/\log F(\mu_\tau;\sigma)$, so that

$$
G(x;\mu_\tau,\sigma) = F(x;\sigma)^{\alpha} = \tau^{\log F(x;\sigma)/\log F(\mu_\tau;\sigma)}
$$

and $G(\mu_\tau;\mu_\tau,\sigma)=\tau$. In
`QEAkash(mu.link = "log", sigma.link = "log", tau = 0.50)`, `tau` is a
fixed scalar in (0, 1). The modeled parameter `mu` is the quantile at
probability level `tau`, and `sigma` is the Akash parameter. Setting
`tau = 0.50` models the conditional median. Under the default links,
`log(mu) = X1 %*% beta1` and `log(sigma) = X2 %*% beta2`. The fitted values
are on the original response scale, and the regression coefficients are on
the log scale.

The family uses the closed-form score and second derivatives of the
log-density. Starting values come from the profile maximum likelihood fit
of the response without covariates (see below), which is exact for
intercept-only models.

## Example 1. Fit a QEAkash distribution

Generate 400 observations with a true 25th percentile of 1.5 and
`sigma = 2`, and estimate both parameters with an intercept-only model.

```r
set.seed(2026)
tau_distribution <- 0.25
distribution_data <- data.frame(
  y = rQEAkash(400, mu = 1.5, sigma = 2, tau = tau_distribution)
)

distribution_fit <- gamlss(
  y ~ 1, sigma.formula = ~ 1,
  family = QEAkash(tau = tau_distribution),
  data = distribution_data,
  control = gamlss.control(trace = FALSE)
)
stopifnot(isTRUE(distribution_fit$converged))
fitted(distribution_fit, what = "mu")[1]
fitted(distribution_fit, what = "sigma")[1]
```

The fitted `mu` is the estimated 25th percentile and is common to all
observations. `qeakash_profile_fit(distribution_data$y, tau_distribution)`
returns the same maximum by profile likelihood.

## Example 2. Fit a QR-QEAkash regression

Generate 600 observations whose conditional median depends on a
continuous covariate `x` and a binary indicator `group`, and whose `sigma`
depends on `x`.

```r
set.seed(2027)
tau_regression <- 0.50
n <- 600
regression_data <- data.frame(
  x = runif(n),
  group = rbinom(n, size = 1, prob = 0.5)
)

mu <- with(regression_data, exp(0.4 + 0.6 * x + 0.3 * group))
sigma <- with(regression_data, exp(0.7 - 0.5 * x))
regression_data$y <- rQEAkash(n, mu = mu, sigma = sigma,
                              tau = tau_regression)

regression_fit <- gamlss(
  y ~ x + group, sigma.formula = ~ x,
  family = QEAkash(tau = tau_regression),
  data = regression_data,
  control = gamlss.control(trace = FALSE, n.cyc = 200)
)
stopifnot(isTRUE(regression_fit$converged))
coef(regression_fit, what = "mu")
coef(regression_fit, what = "sigma")
head(fitted(regression_fit, what = "mu"))
```

Holding the other covariate fixed, a one-unit increase in a predictor of
`mu` multiplies the conditional median by `exp(beta)`, a change of
`100 * (exp(beta) - 1)` percent. The coefficients of the `sigma` submodel
act on `sigma` in the same way, and their effect on dispersion and shape
follows from the properties of the distribution. Use
`residuals(regression_fit)` to obtain the normalized quantile residuals
for diagnostic plots.

## Distribution functions

| Function | Purpose |
|:--|:--|
| `dQEAkash(y, mu, sigma, tau, log = FALSE)` | Density for positive observations, or log-density when `log = TRUE`. At `y = 0` it returns the one-sided limit. |
| `pQEAkash(q, mu, sigma, tau, lower.tail = TRUE, log.p = FALSE)` | Cumulative probability or, with `lower.tail = FALSE`, the survival probability. |
| `qQEAkash(p, mu, sigma, tau)` | Quantile at probability `p` in [0, 1], computed numerically. In particular, `qQEAkash(tau, mu, sigma, tau)` returns `mu` up to rounding. |
| `rQEAkash(n, mu, sigma, tau)` | Generate `n` observations by inverse transform sampling, with common parameters or one `mu` and `sigma` per observation. |
| `sQEAkash(q, mu, sigma, tau)`, `hQEAkash(y, mu, sigma, tau)` | Survival and hazard functions. |
| `qeakash_alpha(mu, sigma, tau)` | Power parameter of the exponentiated Akash distribution. |
| `qeakash_profile_fit(y, tau)` | Maximum likelihood estimates of `mu` and `sigma` for a sample without covariates, by the profile log-likelihood in `sigma`. |

These helpers default to `tau = 0.50`. Pass the same `tau` to the data
generator and fitted family when reproducing a simulation. For a sample
without covariates, the profile log-likelihood in `sigma` does not depend
on `tau`: changing `tau` changes the estimated quantile, not the fitted
distribution.

Right-censored QR-QEAkash regressions are fitted by direct maximization of
the censored log-likelihood with `qr_qeakash_fit()`, in
[Simulations_QR_QEAkash](../Simulations_QR_QEAkash/).
