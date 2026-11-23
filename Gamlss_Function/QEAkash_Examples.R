# Reproducible QEAkash examples
# From the directory containing both files, run:
#   Rscript QEAkash_Examples.R
# Install dependencies once if needed:
#   install.packages(c("gamlss", "gamlss.dist"))

library(gamlss)
source("QEAkash_GAMLSS.R")

# Example 1: Fit a QEAkash distribution with a common 25th percentile.
set.seed(2026)
tau_distribution <- 0.25
true_parameters <- c(mu = 1.50, sigma = 2.00)
distribution_data <- data.frame(
  y = rQEAkash(400, mu = true_parameters[["mu"]],
               sigma = true_parameters[["sigma"]], tau = tau_distribution)
)

distribution_fit <- gamlss(
  y ~ 1, sigma.formula = ~ 1, family = QEAkash(tau = tau_distribution),
  data = distribution_data, control = gamlss.control(trace = FALSE)
)
stopifnot(isTRUE(distribution_fit$converged))
cat("Example 1: Common 25th percentile and sigma\n")
print(cbind(
  true = true_parameters,
  gamlss = c(fitted(distribution_fit, what = "mu")[1],
             fitted(distribution_fit, what = "sigma")[1]),
  profile = unlist(qeakash_profile_fit(distribution_data$y,
                                       tau_distribution)[c("mu", "sigma")])
))

# Example 2: Fit a QR-QEAkash model for the conditional median.
# log(mu) depends on x and group, and log(sigma) depends on x.
set.seed(2027)
tau_regression <- 0.50
n <- 600
regression_data <- data.frame(
  x = stats::runif(n),
  group = stats::rbinom(n, size = 1, prob = 0.5)
)
true_beta1 <- c("(Intercept)" = 0.4, x = 0.6, group = 0.3)
true_beta2 <- c("(Intercept)" = 0.7, x = -0.5)
conditional_quantile <- with(
  regression_data, exp(true_beta1[1] + true_beta1[2] * x +
                         true_beta1[3] * group)
)
conditional_sigma <- with(regression_data, exp(true_beta2[1] + true_beta2[2] * x))
regression_data$y <- rQEAkash(n, mu = conditional_quantile,
                              sigma = conditional_sigma, tau = tau_regression)

regression_fit <- gamlss(
  y ~ x + group, sigma.formula = ~ x,
  family = QEAkash(tau = tau_regression), data = regression_data,
  control = gamlss.control(trace = FALSE, n.cyc = 200)
)
stopifnot(isTRUE(regression_fit$converged))
cat("\nExample 2: Coefficients on the log scale\n")
print(rbind(
  cbind(true = true_beta1, estimated = coef(regression_fit, what = "mu")),
  cbind(true = true_beta2, estimated = coef(regression_fit, what = "sigma"))
))

# Fitted conditional medians are on the original response scale.
cat("\nFirst six fitted conditional medians\n")
print(head(fitted(regression_fit, what = "mu")))

# Normalized quantile residuals should be close to a standard normal sample.
cat("\nShapiro-Wilk test of the quantile residuals\n")
print(stats::shapiro.test(residuals(regression_fit)))
