# QEAkash distribution: data generation and maximum likelihood estimation.
# The distribution is parameterized by its quantile mu_tau at level tau and
# by the Akash parameter sigma:
#   alpha = log(tau) / log F(mu_tau; sigma),  G(y) = F(y; sigma)^alpha,
#   F(y; sigma) = 1 - [1 + sigma y (sigma y + 2)/(sigma^2 + 2)] exp(-sigma y).
# The same estimation routine serves the regression study, where
#   log(mu_tau,t) = x1t' beta1 and log(sigma_t) = x2t' beta2,
# with right censoring (delta_t = 0) allowed. Without covariates, X1 and X2
# are columns of ones. Sourcing this file only defines functions.

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

# One replication of the distribution study. Returns the estimates of
# mu_tau and sigma and the coverage indicators of the 95% Wald intervals
# for log(mu_tau) and log(sigma), or NA values when the fit fails.
qe_replicate_distribution <- function(seed, mu, s, tau, n) {
  set.seed(seed)
  y <- qe_generate(runif(n), mu, s, tau)
  one <- matrix(1, n, 1)
  f <- qe_fit(y, rep(1, n), one, one, tau)
  if (is.null(f)) return(rep(NA_real_, 4))
  c(exp(f$th), as.numeric(abs(f$th - log(c(mu, s))) <= qnorm(0.975) * f$se))
}
