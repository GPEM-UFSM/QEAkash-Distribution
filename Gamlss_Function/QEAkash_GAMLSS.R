# Quantile Exponentiated Akash (QEAkash) distribution for GAMLSS
#
# The Akash distribution (Shanker, 2015) with parameter sigma > 0 has CDF
#   F(x; sigma) = 1 - [1 + sigma x (sigma x + 2) / (sigma^2 + 2)] exp(-sigma x).
# The exponentiated Akash distribution raises it to a power alpha > 0. Fixing
# the probability level tau in (0, 1) and writing mu = mu_tau for the quantile
# of order tau gives
#   alpha = alpha(mu, sigma) = log(tau) / log F(mu; sigma),
#   G(x; mu, sigma) = F(x; sigma)^alpha = tau^(log F(x; sigma) / log F(mu; sigma)),
# so that G(mu; mu, sigma) = tau. Both mu and sigma are positive, and the
# response is positive. See QEAkash_Examples.R for an intercept-only fit and a
# QR-QEAkash regression with submodels for mu and sigma.

# ---------------------------------------------------------------------------
# Numerically stable pieces of the Akash base distribution
# ---------------------------------------------------------------------------

# log(exp(a) + exp(b)), with log(0 + 0) = -Inf.
.qeakash_lse <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

# log F(x; sigma). The Akash distribution is the mixture
#   sigma^2/(sigma^2 + 2) Gamma(1, sigma) + 2/(sigma^2 + 2) Gamma(3, sigma),
# which gives log F accurately near zero. Near the upper tail, log F is
# computed as log(1 - S) from the survival function S to avoid 1 - F
# cancellation.
qeakash_log_akash_cdf <- function(x, sigma) {
  sx <- sigma * x
  lw1 <- 2 * log(sigma) - log(sigma^2 + 2)
  lw2 <- log(2) - log(sigma^2 + 2)
  lF <- .qeakash_lse(lw1 + stats::pgamma(sx, 1, log.p = TRUE),
                     lw2 + stats::pgamma(sx, 3, log.p = TRUE))
  lS <- .qeakash_lse(lw1 - sx,
                     lw2 + stats::pgamma(sx, 3, lower.tail = FALSE, log.p = TRUE))
  ifelse(lS < log(0.5), log1p(-exp(lS)), lF)
}

# log of the Akash density f(x; sigma) for x > 0.
qeakash_log_akash_pdf <- function(x, sigma) {
  3 * log(sigma) - log(sigma^2 + 2) + log1p(x^2) - sigma * x
}

qeakash_validate_tau <- function(tau) {
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) ||
      tau <= 0 || tau >= 1) {
    stop("tau must be a single number strictly between 0 and 1.")
  }
}

qeakash_validate_par <- function(mu, sigma) {
  if (!is.numeric(mu) || !length(mu) || any(!is.finite(mu) | mu <= 0)) {
    stop("mu must contain finite positive quantiles.")
  }
  if (!is.numeric(sigma) || !length(sigma) ||
      any(!is.finite(sigma) | sigma <= 0)) {
    stop("sigma must contain finite positive values.")
  }
}

# alpha(mu, sigma) = log(tau) / log F(mu; sigma), always positive.
qeakash_alpha <- function(mu, sigma, tau = 0.50) {
  qeakash_validate_tau(tau)
  qeakash_validate_par(mu, sigma)
  log(tau) / qeakash_log_akash_cdf(mu, sigma)
}

# ---------------------------------------------------------------------------
# d, p, q, r, survival and hazard functions
# ---------------------------------------------------------------------------

# Density. At y = 0 the one-sided limit is returned: Inf when alpha < 1,
# sigma^3/(sigma^2 + 2) when alpha = 1 and 0 when alpha > 1.
dQEAkash <- function(y, mu, sigma, tau = 0.50, log = FALSE) {
  alpha <- qeakash_alpha(mu, sigma, tau)
  n <- max(length(y), length(alpha), length(sigma))
  y <- rep_len(y, n); alpha <- rep_len(alpha, n); sigma <- rep_len(sigma, n)
  out <- rep(-Inf, n)
  out[is.na(y)] <- NA_real_
  k <- which(!is.na(y) & y > 0 & is.finite(y))
  out[k] <- log(alpha[k]) + qeakash_log_akash_pdf(y[k], sigma[k]) +
    (alpha[k] - 1) * qeakash_log_akash_cdf(y[k], sigma[k])
  k0 <- which(!is.na(y) & y == 0)
  out[k0] <- ifelse(alpha[k0] < 1, Inf, ifelse(alpha[k0] > 1, -Inf,
                    3 * log(sigma[k0]) - log(sigma[k0]^2 + 2)))
  if (log) out else exp(out)
}

# CDF G(q) = F(q; sigma)^alpha. With lower.tail = FALSE it returns the
# survival function, computed as 1 - exp(log G) without cancellation.
pQEAkash <- function(q, mu, sigma, tau = 0.50, lower.tail = TRUE,
                     log.p = FALSE) {
  alpha <- qeakash_alpha(mu, sigma, tau)
  n <- max(length(q), length(alpha), length(sigma))
  q <- rep_len(q, n); alpha <- rep_len(alpha, n); sigma <- rep_len(sigma, n)
  lG <- rep(NA_real_, n)
  lG[!is.na(q) & q <= 0] <- -Inf
  k <- which(!is.na(q) & q > 0)
  lG[k] <- alpha[k] * qeakash_log_akash_cdf(q[k], sigma[k])
  out <- if (lower.tail) lG else log(-expm1(lG))
  if (log.p) out else exp(out)
}

# Survival and hazard functions.
sQEAkash <- function(q, mu, sigma, tau = 0.50, log.p = FALSE) {
  pQEAkash(q, mu, sigma, tau, lower.tail = FALSE, log.p = log.p)
}

hQEAkash <- function(y, mu, sigma, tau = 0.50) {
  exp(dQEAkash(y, mu, sigma, tau, log = TRUE) -
        sQEAkash(y, mu, sigma, tau, log.p = TRUE))
}

# Quantile function Q(p) = F^{-1}(p^(1/alpha); sigma). The Akash quantile
# has no closed form, so log F(x; sigma) = log(p)/alpha is solved for log(x)
# with uniroot. In particular, qQEAkash(tau, mu, sigma, tau) returns mu.
qQEAkash <- function(p, mu, sigma, tau = 0.50, lower.tail = TRUE,
                     log.p = FALSE) {
  alpha <- qeakash_alpha(mu, sigma, tau)
  lp <- if (log.p) p else log(p)
  if (anyNA(lp) || any(lp > 0)) {
    stop("p must contain probabilities between 0 and 1.")
  }
  if (!lower.tail) lp <- log(-expm1(lp))
  n <- max(length(lp), length(alpha), length(mu), length(sigma))
  lp <- rep_len(lp, n); alpha <- rep_len(alpha, n)
  mu <- rep_len(mu, n); sigma <- rep_len(sigma, n)
  vapply(seq_len(n), function(i) {
    if (lp[i] == -Inf) return(0)
    if (lp[i] == 0) return(Inf)
    target <- lp[i] / alpha[i]
    f <- function(t) qeakash_log_akash_cdf(exp(t), sigma[i]) - target
    lo <- log(mu[i]) - 1; hi <- log(mu[i]) + 1
    while (f(lo) > 0 && lo > -700) lo <- lo - 4
    while (f(hi) < 0 && hi < 700) hi <- hi + 4
    exp(stats::uniroot(f, c(lo, hi), tol = 1e-12)$root)
  }, numeric(1))
}

# Random generation by inverse transform sampling. Supply either common
# parameters or one value per observation.
rQEAkash <- function(n, mu, sigma, tau = 0.50) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 0 || n != floor(n)) {
    stop("n must be a nonnegative integer.")
  }
  if (!length(mu) %in% c(1L, n) || !length(sigma) %in% c(1L, n)) {
    stop("mu and sigma must be scalars or have one value per observation.")
  }
  if (n == 0) return(numeric())
  qQEAkash(stats::runif(n), mu, sigma, tau)
}

# ---------------------------------------------------------------------------
# Score and second derivatives of log g(y; mu, sigma) for one observation
# ---------------------------------------------------------------------------
# With L = log F(mu; sigma), M = log F(y; sigma), alpha = log(tau)/L,
# a(x) = f(x; sigma)/F(x; sigma) and b(x) = dF(x; sigma)/dsigma / F(x; sigma):
#   dl/dmu    = alpha_mu (1/alpha + M),       alpha_mu    = -alpha a(mu)/L,
#   dl/dsigma = alpha_sigma (1/alpha + M)     alpha_sigma = -alpha b(mu)/L,
#               + 3/sigma - 2 sigma/(sigma^2 + 2) - y + (alpha - 1) b(y).
# The second derivatives follow by the product and chain rules.
qeakash_derivatives <- function(y, mu, sigma, tau) {
  K <- function(x) (sigma^2 + 2) * (1 + x^2) + 2 * sigma * x + 4
  Ks <- function(x) 2 * sigma * (1 + x^2) + 2 * x
  log_Fs <- function(x) {            # log dF(x; sigma)/dsigma, always positive
    2 * log(sigma) + log(x) - sigma * x + log(K(x)) - 2 * log(sigma^2 + 2)
  }
  L <- qeakash_log_akash_cdf(mu, sigma)
  M <- qeakash_log_akash_cdf(y, sigma)
  A <- log(tau) / L
  a_mu <- exp(qeakash_log_akash_pdf(mu, sigma) - L)
  b_mu <- exp(log_Fs(mu) - L)
  b_y <- exp(log_Fs(y) - M)
  ds <- 3 / sigma - 2 * sigma / (sigma^2 + 2)
  L_mm <- a_mu * (2 * mu / (1 + mu^2) - sigma) - a_mu^2
  L_ms <- a_mu * (ds - mu) - a_mu * b_mu
  L_ss <- b_mu * (2 / sigma - mu + Ks(mu) / K(mu) - 4 * sigma / (sigma^2 + 2)) -
    b_mu^2
  M_ss <- b_y * (2 / sigma - y + Ks(y) / K(y) - 4 * sigma / (sigma^2 + 2)) -
    b_y^2
  A_m <- -A * a_mu / L
  A_s <- -A * b_mu / L
  A_mm <- A * (2 * a_mu^2 / L^2 - L_mm / L)
  A_ms <- A * (2 * a_mu * b_mu / L^2 - L_ms / L)
  A_ss <- A * (2 * b_mu^2 / L^2 - L_ss / L)
  V <- 1 / A + M
  list(
    dldm = A_m * V,
    dldd = ds - y + A_s * V + (A - 1) * b_y,
    d2ldm2 = A_mm * V - A_m^2 / A^2,
    d2ldmdd = A_ms * V - A_m * A_s / A^2 + A_m * b_y,
    d2ldd2 = -3 / sigma^2 + 2 * (sigma^2 - 2) / (sigma^2 + 2)^2 + A_ss * V -
      A_s^2 / A^2 + 2 * A_s * b_y + (A - 1) * M_ss
  )
}

# ---------------------------------------------------------------------------
# Maximum likelihood for a sample without covariates (profile likelihood)
# ---------------------------------------------------------------------------
# For fixed sigma, the score in alpha vanishes at alpha(sigma) = -n/M(sigma),
# with M(sigma) = sum log F(y_t; sigma). The corresponding quantile is
#   mu(sigma) = F^{-1}(tau^(-M(sigma)/n); sigma),
# and the profile log-likelihood
#   l_p(sigma) = n log(-n/M) + 3n log(sigma) - n log(sigma^2 + 2)
#                + sum log(1 + y^2) - sigma sum y - n - M
# does not depend on tau. sigma is found by a one-dimensional search.
qeakash_profile_loglik <- function(sigma, y) {
  n <- length(y)
  vapply(sigma, function(s) {
    M <- sum(qeakash_log_akash_cdf(y, s))
    if (!is.finite(M) || M >= 0) return(-Inf)
    n * log(-n / M) + 3 * n * log(s) - n * log(s^2 + 2) + sum(log1p(y^2)) -
      s * sum(y) - n - M
  }, numeric(1))
}

qeakash_profile_fit <- function(y, tau = 0.50) {
  qeakash_validate_tau(tau)
  if (!is.numeric(y) || length(y) < 2L || any(!is.finite(y) | y <= 0)) {
    stop("y must contain at least two finite positive observations.")
  }
  f <- function(ls) -qeakash_profile_loglik(exp(ls), y)
  # Akash scale: the mean is close to 1/sigma for large sigma and to
  # 3/sigma for small sigma, so the search grid is centred on log(1/mean(y)).
  grid <- log(1 / mean(y)) + seq(-8, 8, by = 0.5)
  vals <- vapply(grid, function(g) tryCatch(f(g), error = function(e) Inf),
                 numeric(1))
  k <- which.min(vals)
  opt <- stats::optimize(f, grid[k] + c(-0.5, 0.5), tol = 1e-12)
  sigma <- exp(opt$minimum)
  M <- sum(qeakash_log_akash_cdf(y, sigma))
  alpha <- -length(y) / M
  # mu(sigma) = F^{-1}(tau^(1/alpha); sigma), solved on the log scale
  target <- log(tau) / alpha
  g <- function(t) qeakash_log_akash_cdf(exp(t), sigma) - target
  lo <- log(stats::median(y)) - 1; hi <- log(stats::median(y)) + 1
  while (g(lo) > 0 && lo > -700) lo <- lo - 4
  while (g(hi) < 0 && hi < 700) hi <- hi + 4
  mu <- exp(stats::uniroot(g, c(lo, hi), tol = 1e-12)$root)
  list(mu = mu, sigma = sigma, alpha = alpha, loglik = -opt$objective)
}

# ---------------------------------------------------------------------------
# GAMLSS family
# ---------------------------------------------------------------------------
# Construct a family with two modeled parameters, mu (the quantile of order
# tau) and sigma. The quantile level is captured by the family and stays
# fixed throughout fitting. Both default links are logarithmic, as in the
# QR-QEAkash model: log(mu) = X1 beta1 and log(sigma) = X2 beta2.
QEAkash <- function(mu.link = "log", sigma.link = "log", tau = 0.50) {
  qeakash_validate_tau(tau)
  force(tau)
  mstats <- gamlss.dist::checklink(
    "mu.link", "QEAkash", substitute(mu.link),
    c("log", "identity", "inverse", "own")
  )
  dstats <- gamlss.dist::checklink(
    "sigma.link", "QEAkash", substitute(sigma.link),
    c("log", "identity", "inverse", "own")
  )

  # GAMLSS evaluates some family components in its own fitting environment.
  # Each evaluator is embedded as a function object in the callback body so
  # that the fixed quantile level remains available after that extraction.
  bind_callback <- function(evaluator) {
    callback <- function(y, mu, sigma, ...) NULL
    body(callback) <- as.call(list(evaluator, quote(y), quote(mu),
                                   quote(sigma)))
    callback
  }
  derivative <- function(name, floor = FALSE) {
    bind_callback(function(y, mu, sigma) {
      value <- qeakash_derivatives(y, mu, sigma, tau)[[name]]
      # Negative curvature floor for the RS/CG working weights.
      if (floor) value <- ifelse(value < -1e-15, value, -1e-15)
      value[!is.finite(value)] <- if (floor) -1e-15 else 0
      value
    })
  }
  deviance <- bind_callback(function(y, mu, sigma) {
    -2 * dQEAkash(y, mu, sigma, tau, log = TRUE)
  })
  quantile_residuals <- bind_callback(function(y, mu, sigma) {
    stats::qnorm(pQEAkash(y, mu, sigma, tau))
  })
  initial <- function(y) {
    fit <- tryCatch(qeakash_profile_fit(y, tau), error = function(e) NULL)
    if (is.null(fit)) {
      list(mu = stats::quantile(y, tau, names = FALSE), sigma = 1 / mean(y))
    } else {
      fit
    }
  }
  mu_initial <- function(y) rep(initial(y)$mu, length(y))
  sigma_initial <- function(y) rep(initial(y)$sigma, length(y))

  structure(
    list(
      family = c("QEAkash", "Quantile Exponentiated Akash"),
      parameters = list(mu = TRUE, sigma = TRUE),
      nopar = 2,
      type = "Continuous",
      tau = tau,
      mu.link = as.character(substitute(mu.link)),
      sigma.link = as.character(substitute(sigma.link)),
      mu.linkfun = mstats$linkfun,
      sigma.linkfun = dstats$linkfun,
      mu.linkinv = mstats$linkinv,
      sigma.linkinv = dstats$linkinv,
      mu.dr = mstats$mu.eta,
      sigma.dr = dstats$mu.eta,
      dldm = derivative("dldm"),
      d2ldm2 = derivative("d2ldm2", floor = TRUE),
      dldd = derivative("dldd"),
      d2ldd2 = derivative("d2ldd2", floor = TRUE),
      d2ldmdd = derivative("d2ldmdd"),
      G.dev.incr = deviance,
      rqres = as.expression(body(quantile_residuals)),
      mu.initial = as.expression(bquote(mu <- .(mu_initial)(y))),
      sigma.initial = as.expression(bquote(sigma <- .(sigma_initial)(y))),
      mu.valid = function(mu) all(is.finite(mu) & mu > 0),
      sigma.valid = function(sigma) all(is.finite(sigma) & sigma > 0),
      y.valid = function(y) all(is.finite(y) & y > 0)
    ),
    class = c("gamlss.family", "family")
  )
}
