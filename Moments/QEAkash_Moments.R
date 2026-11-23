# QEAkash raw moments
#
# This script computes E[X^r] for nonnegative integer orders r under the
# Quantile Exponentiated Akash (QEAkash) distribution, by adaptive
# quadrature and by the binomial series of the density. The series uses the
# Rmpfr package for multiple-precision arithmetic when it is installed.
#
# For a single moment, source this file in a separate environment:
#   qeakash <- new.env()
#   sys.source("QEAkash_Moments.R", envir = qeakash)
#   qeakash$qeakash_moment_integral(r = 1, mu_tau = 2.5, sigma = 1.8, tau = 0.5)
#
# Run Rscript QEAkash_Moments.R to reproduce the full comparison. The script
# prints the results and writes qeakash_moment_comparison.csv and
# qeakash_moment_comparison.tex to the current working directory. The LaTeX
# table uses booktabs commands.
#
# With alpha = log(tau)/log F(mu_tau; sigma), W(x) = 1 + A x + B x^2,
# A = 2 sigma/(sigma^2 + 2), B = sigma^2/(sigma^2 + 2) and c_j = (j + 1) sigma,
#   E[X^r] = alpha sigma^3/(sigma^2 + 2) sum_{j >= 0} (-1)^j C(alpha - 1, j)
#            sum_{m = 0}^{2j} b_{j,m} Lambda_r(j, m),
#   Lambda_r(j, m) = (r + m)!/c_j^(r + m + 1) + (r + m + 2)!/c_j^(r + m + 3),
# where b_{j,m} is the coefficient of x^m in W(x)^j. The series converges
# absolutely for every alpha > 0. It is finite only when alpha - 1 is a
# nonnegative integer. The comparison uses direct quadrature as its
# numerical reference.

# Akash distribution and density functions.
qeakash_akash_cdf <- function(x, sigma) {
  u <- sigma * x * (sigma * x + 2) / (sigma^2 + 2)
  ue <- ifelse(sigma * x > 745, 0, u * exp(-sigma * x))
  # 1 - (1 + u) e^(-sigma x) without cancellation for small x.
  -expm1(-sigma * x) - ue
}

qeakash_akash_pdf <- function(x, sigma) {
  sigma^3 / (sigma^2 + 2) * (1 + x^2) * exp(-sigma * x)
}

# Convert the quantile mu_tau and its probability level tau to the positive
# power parameter alpha. All inputs must be scalars.
qeakash_alpha <- function(mu_tau, sigma, tau) {
  stopifnot(length(mu_tau) == 1L, length(sigma) == 1L, length(tau) == 1L,
            is.finite(mu_tau), is.finite(sigma), is.finite(tau),
            mu_tau > 0, sigma > 0, tau > 0, tau < 1)
  log(tau) / log(qeakash_akash_cdf(mu_tau, sigma))
}

# QEAkash density for given alpha and sigma.
qeakash_pdf <- function(x, alpha, sigma) {
  alpha * qeakash_akash_pdf(x, sigma) *
    qeakash_akash_cdf(x, sigma)^(alpha - 1)
}

# QEAkash quantile, used to split the integration range. It solves
# F(x; sigma) = F(mu_tau; sigma)^(log(p)/log(tau)).
qeakash_quantile <- function(p, mu_tau, sigma, tau) {
  target <- qeakash_akash_cdf(mu_tau, sigma)^(log(p) / log(tau))
  upper <- mu_tau + 1
  while (qeakash_akash_cdf(upper, sigma) < target) upper <- upper * 2
  stats::uniroot(function(x) qeakash_akash_cdf(x, sigma) - target,
                 lower = 0, upper = upper, tol = 1e-12)$root
}

# Adaptive Gauss-Kronrod quadrature of x^r g(x) over (0, Inf). The range is
# split at the median. When alpha < 1 the density behaves as x^(alpha - 1)
# near zero, and the substitution v = x^alpha removes that singularity.
qeakash_moment_integral <- function(r, mu_tau, sigma, tau, rel.tol = 1e-11) {
  stopifnot(r >= 0)
  alpha <- qeakash_alpha(mu_tau, sigma, tau)
  f <- function(x) x^r * qeakash_pdf(x, alpha, sigma)
  cut <- qeakash_quantile(0.5, mu_tau, sigma, tau)
  if (alpha >= 1) {
    lower_part <- stats::integrate(f, 0, cut, rel.tol = rel.tol,
                                   subdivisions = 5000L)$value
  } else {
    lower_part <- stats::integrate(function(v) {
      value <- f(v^(1 / alpha)) * (1 / alpha) * v^(1 / alpha - 1)
      value[!is.finite(value)] <- 0
      value
    }, 0, cut^alpha, rel.tol = rel.tol, subdivisions = 5000L)$value
  }
  lower_part + stats::integrate(f, cut, Inf, rel.tol = rel.tol,
                                subdivisions = 5000L)$value
}

# Binomial series truncated at j = J, for one or more orders r. The
# coefficients of W(x)^j follow the convolution recurrence
#   b_{0,0} = 1, b_{j,m} = b_{j-1,m} + A b_{j-1,m-1} + B b_{j-1,m-2}.
# For j < alpha the weights (-1)^j C(alpha - 1, j) alternate in sign and
# grow, so the sum is accumulated with precision_bits bits when Rmpfr is
# available. Set precision_bits = NULL to use double precision.
qeakash_moment_series <- function(r, mu_tau, sigma, tau, J = 200L,
                                  precision_bits = 512L) {
  stopifnot(all(r >= 0), J >= 0)
  use_mpfr <- !is.null(precision_bits) &&
    requireNamespace("Rmpfr", quietly = TRUE)
  if (!use_mpfr) {
    if (!is.null(precision_bits)) {
      warning("Rmpfr is not installed; using double precision.")
    }
    return(qeakash_moment_series_double(r, mu_tau, sigma, tau, J))
  }
  num <- function(x) Rmpfr::mpfr(x, precision_bits)
  zeros <- function(n) Rmpfr::mpfr(rep(0, n), precision_bits)
  fact <- function(k) factorial(Rmpfr::mpfr(k, precision_bits))
  s <- num(sigma)
  m0 <- num(mu_tau)
  tt <- num(tau)
  F_mu <- 1 - (1 + s * m0 * (s * m0 + 2) / (s^2 + 2)) * exp(-s * m0)
  alpha <- log(tt) / log(F_mu)
  A <- 2 * s / (s^2 + 2)
  B <- s^2 / (s^2 + 2)
  constant <- alpha * s^3 / (s^2 + 2)
  factorials <- fact(0:(2 * J + max(r) + 2))   # factorials[k + 1] = k!
  b <- num(1)
  weight <- num(1)
  total <- zeros(length(r))
  for (j in 0:J) {
    if (j > 0) {
      weight <- weight * (j - alpha) / j
      n0 <- length(b)
      z <- zeros(n0 + 2)
      z[1:n0] <- z[1:n0] + b
      z[2:(n0 + 1)] <- z[2:(n0 + 1)] + A * b
      z[3:(n0 + 2)] <- z[3:(n0 + 2)] + B * b
      b <- z
    }
    m <- 0:(length(b) - 1L)
    c_j <- (j + 1) * s
    for (k in seq_along(r)) {
      lambda <- factorials[r[k] + m + 1] / c_j^(r[k] + m + 1) +
        factorials[r[k] + m + 3] / c_j^(r[k] + m + 3)
      total[k] <- total[k] + weight * sum(b * lambda)
    }
  }
  as.numeric(constant * total)
}

# Double-precision version of the same truncated series. The coefficients
# b_{j,m} are positive and reach values far below the smallest double, and
# the factorials overflow beyond 170!, so both are kept on the log scale.
qeakash_moment_series_double <- function(r, mu_tau, sigma, tau, J = 200L) {
  lse3 <- function(a, b, c) {
    m <- pmax(a, b, c)
    out <- m + log(exp(a - m) + exp(b - m) + exp(c - m))
    out[!is.finite(m)] <- -Inf
    out
  }
  alpha <- qeakash_alpha(mu_tau, sigma, tau)
  log_A <- log(2 * sigma / (sigma^2 + 2))
  log_B <- log(sigma^2 / (sigma^2 + 2))
  log_b <- 0
  weight <- 1
  total <- numeric(length(r))
  for (j in 0:J) {
    if (j > 0) {
      weight <- weight * (j - alpha) / j
      log_b <- lse3(c(log_b, -Inf, -Inf), log_A + c(-Inf, log_b, -Inf),
                    log_B + c(-Inf, -Inf, log_b))
    }
    m <- seq_along(log_b) - 1L
    log_c <- log((j + 1) * sigma)
    for (k in seq_along(r)) {
      total[k] <- total[k] + weight * sum(
        exp(log_b + lfactorial(r[k] + m) - (r[k] + m + 1) * log_c) +
          exp(log_b + lfactorial(r[k] + m + 2) - (r[k] + m + 3) * log_c))
    }
  }
  alpha * sigma^3 / (sigma^2 + 2) * total
}

# Closed-form Akash moments, reached when alpha = 1, that is,
# tau = F(mu_tau; sigma).
qeakash_moment_akash <- function(r, sigma) {
  factorial(r) * (sigma^2 + (r + 1) * (r + 2)) / (sigma^r * (sigma^2 + 2))
}

# Compare quadrature and the truncated series for each case and order.
# cases is a data frame with columns mu_tau and sigma. D_r is the absolute
# difference between the two approximations.
qeakash_compare_moments <- function(cases, taus = c(0.10, 0.25, 0.50, 0.75, 0.90),
                                    moments = 1:4, J = 200L,
                                    precision_bits = 512L) {
  stopifnot(all(c("mu_tau", "sigma") %in% names(cases)))
  rows <- list()
  for (i in seq_len(nrow(cases))) {
    for (tau in taus) {
      mu_tau <- cases$mu_tau[i]
      sigma <- cases$sigma[i]
      series <- qeakash_moment_series(moments, mu_tau, sigma, tau, J,
                                      precision_bits)
      integral <- vapply(moments, qeakash_moment_integral, numeric(1),
                         mu_tau = mu_tau, sigma = sigma, tau = tau)
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = i, mu_tau = mu_tau, sigma = sigma, tau = tau,
        alpha = qeakash_alpha(mu_tau, sigma, tau), r = moments,
        moment_integral = integral, moment_series = series,
        abs_difference = abs(integral - series)
      )
    }
  }
  do.call(rbind, rows)
}

# LaTeX rows with four decimals for the moments and D_r in scientific
# notation, two blocks per scenario (r = 1, 2 and r = 3, 4).
qeakash_latex_rows <- function(comparison) {
  fmt_moment <- function(x) formatC(x, format = "f", digits = 4)
  fmt_difference <- function(x) {
    if (x == 0) return("$0$")
    e <- floor(log10(x)); m <- x / 10^e
    if (round(m, 2) >= 10) { m <- m / 10; e <- e + 1 }
    sprintf("$%s\\times10^{%d}$", formatC(m, format = "f", digits = 2), e)
  }
  out <- character()
  for (i in unique(comparison$scenario)) {
    ci <- comparison[comparison$scenario == i, ]
    out <- c(out, "\\midrule", sprintf(
      "\\multicolumn{7}{l}{Scenario %d: $\\mu_\\tau = %s$, $\\sigma = %s$} \\\\",
      i, format(ci$mu_tau[1]), format(ci$sigma[1])))
    for (block in list(1:2, 3:4)) {
      if (!all(block %in% ci$r)) next
      out <- c(out, sprintf("$\\tau$ & %s \\\\", paste(sprintf(
        "$\\mu_%d'^{\\,I}$ & $\\mu_%d'^{\\,S}$ & $D_%d$", block, block, block),
        collapse = " & ")))
      for (tau in unique(ci$tau)) {
        d <- ci[ci$tau == tau, ]
        cells <- unlist(lapply(block, function(r) {
          x <- d[d$r == r, ]
          c(sprintf("$%s$", fmt_moment(x$moment_integral)),
            sprintf("$%s$", fmt_moment(x$moment_series)),
            fmt_difference(x$abs_difference))
        }))
        out <- c(out, sprintf("$%s$ & %s \\\\", formatC(tau, format = "f",
                                                         digits = 2),
                              paste(cells, collapse = " & ")))
      }
    }
  }
  out
}

if (identical(environment(), globalenv())) {
  cases <- data.frame(mu_tau = c(0.5, 1.5, 2.5, 0.5, 1.5, 2.5, 0.5, 1.5),
                      sigma = c(0.9, 0.9, 0.9, 1.8, 1.8, 1.8, 2.9, 2.9))
  if (!requireNamespace("Rmpfr", quietly = TRUE)) {
    message("Rmpfr is not installed; the series uses double precision. ",
            "Install it with install.packages(\"Rmpfr\") to reproduce the ",
            "archived comparison.")
  }
  comparison <- qeakash_compare_moments(cases, moments = 1:4, J = 200L)
  options(width = 160)
  print(format(comparison, digits = 6), row.names = FALSE)
  cat(sprintf("\nRange of D_r: %.2e to %.2e\n", min(comparison$abs_difference),
              max(comparison$abs_difference)))
  utils::write.csv(comparison, "qeakash_moment_comparison.csv",
                   row.names = FALSE)
  writeLines(c("\\begin{tabular}{c rrr rrr}", "\\toprule",
               qeakash_latex_rows(comparison), "\\bottomrule",
               "\\end{tabular}"),
             "qeakash_moment_comparison.tex")
  cat("Wrote qeakash_moment_comparison.csv and qeakash_moment_comparison.tex\n")
}
