# QR-QEAkash regression: data generation, right censoring and maximum
# likelihood estimation. For a fixed quantile level tau,
#   log(mu_tau,t) = x1t' beta1,  log(sigma_t) = x2t' beta2,
# and the response follows the QEAkash distribution
#   G(y) = F(y; sigma)^alpha,  alpha = log(tau) / log F(mu_tau; sigma),
#   F(y; sigma) = 1 - [1 + sigma y (sigma y + 2)/(sigma^2 + 2)] exp(-sigma y).
# With right censoring, failures (d = 1) contribute log g(y) and censored
# observations (d = 0) contribute log(1 - G(y)) to the log-likelihood.
# Sourcing this file only defines functions.

# log(exp(a) + exp(b)).
qe_lse <- function(a, b) {
  m <- pmax(a, b)
  m + log1p(exp(-abs(a - b)))
}

# log(-log F(x; sigma)) for the Akash distribution. F is computed from the
# mixture sigma^2/(sigma^2+2) Exp(sigma) + 2/(sigma^2+2) Gamma(3, sigma), and
# log(-log F) remains accurate when F is close to zero and close to one.
qe_lmlF <- function(x, s) {
  u <- s * x
  ls2 <- 2 * log(s)
  lc <- log(s^2 + 2)
  lS <- qe_lse(ls2 - u, log(2) + pgamma(u, 3, lower.tail = FALSE,
                                        log.p = TRUE)) - lc
  lF <- qe_lse(ls2 + pgamma(u, 1, log.p = TRUE),
               log(2) + pgamma(u, 3, log.p = TRUE)) - lc
  i <- lS < log(0.5)
  out <- numeric(length(lS))
  out[!i] <- log(-lF[!i])
  out[i] <- log(-log1p(-exp(lS[i])))
  j <- lS < -700
  out[j] <- lS[j]
  out
}

qe_ltau <- function(tau) log(-log(tau))

# Power parameter alpha(mu_tau, sigma).
qe_alpha <- function(mu, s, tau) exp(qe_ltau(tau) - qe_lmlF(mu, s))

# Log-likelihood contributions: log g(y) for failures (d = 1) and
# log(1 - G(y)) for right-censored observations (d = 0).
qe_loglik_obs <- function(y, d, mu, s, tau) {
  lmy <- qe_lmlF(y, s)
  la <- qe_ltau(tau) - qe_lmlF(mu, s)                  # log(alpha)
  A <- -exp(la + lmy)                                  # log G(y) = alpha log F(y)
  lf <- 3 * log(s) - log(s^2 + 2) + log1p(y^2) - s * y # log f(y; sigma)
  ifelse(d == 1, la + lf + A + exp(lmy), log(-expm1(A)))
}

# Inverse transform sampling on the log scale:
#   log(-log F(T; sigma)) = log(-log U) - log(alpha),
# solved by bisection in log(T).
qe_generate <- function(u, mu, s, tau) {
  n <- length(u)
  s <- rep_len(s, n)
  mu <- rep_len(mu, n)
  target <- log(-log(u)) - (qe_ltau(tau) - qe_lmlF(mu, s))
  lo <- rep(-200, n)
  hi <- log(mu) + 1
  repeat {
    k <- qe_lmlF(exp(lo), s) <= target
    if (!any(k)) break
    lo[k] <- pmax(lo[k] * 1.5, -740)
  }
  repeat {
    k <- qe_lmlF(exp(hi), s) > target
    if (!any(k)) break
    hi[k] <- hi[k] + 2
  }
  for (it in 1:90) {
    md <- (lo + hi) / 2
    k <- qe_lmlF(exp(md), s) > target
    lo[k] <- md[k]
    hi[!k] <- md[!k]
  }
  exp((lo + hi) / 2)
}

# Derivatives of L = log(-log F(x; s)) in x and in s, on the log scale.
# dF/ds = s^2 x (s^2 + 6 + (s^2 + 2) x^2 + 2 s x) exp(-s x)/(s^2 + 2)^2 > 0.
qe_dL <- function(x, s, L) {
  lc <- log(s^2 + 2)
  lF <- -exp(L)
  lf <- 3 * log(s) - lc + log1p(x^2) - s * x
  lFs <- 2 * log(s) + log(x) - s * x +
    log(s^2 + 6 + (s^2 + 2) * x^2 + 2 * s * x) - 2 * lc
  list(x = -exp(lf - lF - L), s = -exp(lFs - lF - L))
}

# Score of each observation with respect to eta1 = log(mu) and
# eta2 = log(sigma), with and without censoring.
qe_score_obs <- function(y, d, mu, s, tau) {
  Ly <- qe_lmlF(y, s)
  Lm <- qe_lmlF(mu, s)
  lA <- qe_ltau(tau) - Lm + Ly
  mA <- exp(lA)
  A <- -mA                                             # log G(y)
  dy <- qe_dL(y, s, Ly)
  dm <- qe_dL(mu, s, Lm)
  am <- -dm$x * mu                                     # dlog(alpha)/deta1
  # Weight of censored observations, -A G/(1 - G), computed so that it stays
  # finite when G(y) falls below the smallest double.
  w <- exp(lA - ifelse(mA > 30, mA + log1p(-exp(-mA)), log(expm1(mA))))
  u1 <- ifelse(d == 1, am * (1 + A), w * am)
  us <- ifelse(d == 1,
               -dm$s + 3 / s - 2 * s / (s^2 + 2) - y + A * (dy$s - dm$s) +
                 exp(Ly) * dy$s,
               w * (dy$s - dm$s))
  cbind(u1, us * s)
}

# Negative log-likelihood and its gradient in theta = (beta1, beta2).
qe_nll <- function(th, y, d, X1, X2, tau) {
  p1 <- ncol(X1)
  mu <- exp(drop(X1 %*% th[1:p1]))
  s <- exp(drop(X2 %*% th[-(1:p1)]))
  v <- qe_loglik_obs(y, d, mu, s, tau)
  if (!all(is.finite(v))) return(1e100)
  -sum(v)
}

qe_gr_nll <- function(th, y, d, X1, X2, tau) {
  p1 <- ncol(X1)
  mu <- exp(drop(X1 %*% th[1:p1]))
  s <- exp(drop(X2 %*% th[-(1:p1)]))
  u <- qe_score_obs(y, d, mu, s, tau)
  -c(crossprod(X1, u[, 1]), crossprod(X2, u[, 2]))
}

# Covariates are standardized internally (Z = X M, beta = M gamma); the
# maximum likelihood estimate is invariant to this reparameterization.
qe_standardize <- function(X) {
  p <- ncol(X)
  M <- diag(p)
  if (p > 1) {
    m <- colMeans(X[, -1, drop = FALSE])
    dp <- apply(X[, -1, drop = FALSE], 2, sd)
    if (any(!is.finite(dp) | dp == 0)) return(NULL)
    M[1, -1] <- -m / dp
    M[cbind(2:p, 2:p)] <- 1 / dp
  }
  list(Z = X %*% M, M = M)
}

# Maximum likelihood fit. Starting values: least squares of log(y) for beta1
# (the sample quantile of order tau without covariates) and a
# one-dimensional profile search for the intercept of log(sigma). nlminb
# with the analytic score is followed by up to three Newton steps with the
# observed information (numerical Jacobian of the analytic score). The
# maximum is accepted when the Newton decrement g' H^{-1} g is below 1e-6
# and the observed information is positive definite. Returns NULL on
# failure, otherwise the estimates and their standard errors, from the
# inverse observed information, on the scale of the linear predictors.
qe_fit <- function(y, d, X1, X2, tau) {
  P1 <- qe_standardize(X1)
  P2 <- qe_standardize(X2)
  if (is.null(P1) || is.null(P2)) return(NULL)
  Z1 <- P1$Z; Z2 <- P2$Z
  p1 <- ncol(Z1); p2 <- ncol(Z2)
  f <- function(t) qe_nll(t, y, d, Z1, Z2, tau)
  g <- function(t) {
    v <- qe_gr_nll(t, y, d, Z1, Z2, tau)
    if (all(is.finite(v))) v else rep(1e10, length(v))
  }
  g1 <- if (p1 == 1) log(quantile(y, tau, names = FALSE)) else
    coef(lm.fit(Z1, log(y)))
  profile <- function(l) f(c(g1, l, rep(0, p2 - 1)))
  grid <- seq(-12, 6, by = 1)
  l0 <- grid[which.min(vapply(grid, profile, numeric(1)))]
  l0 <- optimize(profile, l0 + c(-1, 1))$minimum
  ctl <- list(iter.max = 1000, eval.max = 2000)
  o <- tryCatch(nlminb(c(g1, l0, rep(0, p2 - 1)), f, g, control = ctl),
                error = function(e) NULL)
  for (k in 1:2) if (!is.null(o) && o$convergence != 0)
    o <- tryCatch(nlminb(o$par, f, g, control = ctl), error = function(e) NULL)
  if (is.null(o) || o$convergence != 0 || o$objective >= 1e99) return(NULL)
  th <- o$par
  for (k in 1:4) {
    H <- tryCatch(numDeriv::jacobian(g, th), error = function(e) NULL)
    if (is.null(H) || any(!is.finite(H))) return(NULL)
    H <- (H + t(H)) / 2
    if (min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) <= 0)
      return(NULL)
    gg <- g(th)
    step <- solve(H, gg)
    dec <- sum(gg * step)
    if (dec < 1e-10 || k == 4) break
    new <- th - step
    if (f(new) <= f(th)) th <- new else break
  }
  if (!is.finite(dec) || dec > 1e-6) return(NULL)
  M <- matrix(0, p1 + p2, p1 + p2)
  M[1:p1, 1:p1] <- P1$M
  M[-(1:p1), -(1:p1)] <- P2$M
  V <- M %*% solve(H) %*% t(M)
  list(th = drop(M %*% th), se = sqrt(diag(V)), vcov = V,
       loglik = -qe_nll(drop(M %*% th), y, d, X1, X2, tau))
}

# Maximum likelihood fit of a QR-QEAkash model, with or without right
# censoring. y: positive responses; delta: 1 for observed failures and 0
# for right-censored observations; X1, X2: design matrices of log(mu_tau)
# and log(sigma), including a column of ones for intercepts. Returns NULL
# when the fit fails, otherwise the coefficients, their standard errors from
# the inverse observed information, the covariance matrix and the maximized
# log-likelihood.
qr_qeakash_fit <- function(y, delta = rep(1, length(y)), X1, X2, tau = 0.5) {
  X1 <- as.matrix(X1); X2 <- as.matrix(X2)
  stopifnot(length(y) == nrow(X1), nrow(X1) == nrow(X2),
            length(delta) == length(y), all(delta %in% c(0, 1)),
            all(is.finite(y) & y > 0), tau > 0, tau < 1)
  f <- qe_fit(y, delta, X1, X2, tau)
  if (is.null(f)) return(NULL)
  nm1 <- colnames(X1); nm2 <- colnames(X2)
  if (is.null(nm1)) nm1 <- paste0("x", seq_len(ncol(X1)))
  if (is.null(nm2)) nm2 <- paste0("x", seq_len(ncol(X2)))
  nm <- c(paste0("mu:", nm1), paste0("sigma:", nm2))
  dimnames(f$vcov) <- list(nm, nm)
  list(coefficients = setNames(f$th, nm), se = setNames(f$se, nm),
       vcov = f$vcov, loglik = f$loglik, n = length(y),
       censored = sum(delta == 0), tau = tau)
}

# One replication of the regression study. The covariates of n units are
# drawn with replacement from the rows of the design pool, failure times
# are generated from the QEAkash model and censoring times from U(0, nu)
# (no censoring when nu = Inf). Returns the estimates, the coverage
# indicators of the 95% Wald intervals and the censoring proportion, or NA
# values when the fit fails.
qe_replicate_regression <- function(seed, scenario, nu, n, tau = 0.5) {
  set.seed(seed)
  i <- sample.int(nrow(scenario$X1), n, replace = TRUE)
  X1 <- scenario$X1[i, , drop = FALSE]
  X2 <- scenario$X2[i, , drop = FALSE]
  Tt <- qe_generate(runif(n), exp(drop(X1 %*% scenario$beta1)),
                    exp(drop(X2 %*% scenario$beta2)), tau)
  C <- if (is.finite(nu)) runif(n, 0, nu) else rep(Inf, n)
  y <- pmin(Tt, C)
  d <- as.integer(Tt <= C)
  f <- qe_fit(y, d, X1, X2, tau)
  p <- length(scenario$beta1) + length(scenario$beta2)
  if (is.null(f)) return(rep(NA_real_, 2 * p + 1))
  c(f$th, as.numeric(abs(f$th - c(scenario$beta1, scenario$beta2)) <=
                       qnorm(0.975) * f$se), 1 - mean(d))
}

# Upper limit nu of C ~ U(0, nu) for a target censoring proportion. Since
# P(C < T) = E[min(T, nu)]/nu, nu solves mean(min(T, nu))/nu = target over
# m failure times generated with a fixed seed.
qe_calibrate_nu <- function(scenario, target, m = 2e5, seed = 7, tau = 0.5) {
  if (target == 0) return(Inf)
  set.seed(seed)
  i <- sample.int(nrow(scenario$X1), m, replace = TRUE)
  Tt <- qe_generate(runif(m), exp(drop(scenario$X1[i, ] %*% scenario$beta1)),
                    exp(drop(scenario$X2[i, ] %*% scenario$beta2)), tau)
  uniroot(function(nu) mean(pmin(Tt, nu)) / nu - target, c(1, 1e7),
          tol = 1e-8)$root
}
