library(shiny)
library(bslib)
library(ggplot2)
library(plotly)

# ============================================================================
# QEAkash (Quantile Exponentiated Akash) DISTRIBUTION FUNCTIONS
# ============================================================================
# Akash distribution with parameter sigma > 0:
#   F(x; sigma) = 1 - [1 + sigma x (sigma x + 2)/(sigma^2 + 2)] exp(-sigma x),
#   f(x; sigma) = sigma^3/(sigma^2 + 2) (1 + x^2) exp(-sigma x).
# For a fixed probability level tau in (0, 1) and the quantile mu = mu_tau,
#   alpha = log(tau) / log F(mu; sigma) > 0,
#   G(x; mu, sigma) = F(x; sigma)^alpha,   G(mu; mu, sigma) = tau,
#   g(x; mu, sigma) = alpha f(x; sigma) F(x; sigma)^(alpha - 1).

# log(exp(a) + exp(b)), with log(0 + 0) = -Inf.
lse <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

# log F(x; sigma) and log S(x; sigma) of the Akash distribution, from the
# mixture sigma^2/(sigma^2 + 2) Exp(sigma) + 2/(sigma^2 + 2) Gamma(3, sigma).
akash_log_cdf_surv <- function(x, sigma) {
  sx <- sigma * x
  lw1 <- 2 * log(sigma) - log(sigma^2 + 2)
  lw2 <- log(2) - log(sigma^2 + 2)
  lF <- lse(lw1 + pgamma(sx, 1, log.p = TRUE), lw2 + pgamma(sx, 3, log.p = TRUE))
  lS <- lse(lw1 - sx, lw2 + pgamma(sx, 3, lower.tail = FALSE, log.p = TRUE))
  list(lF = lF, lS = lS)
}

# log F(x; sigma), accurate near zero and in the upper tail.
akash_log_cdf <- function(x, sigma) {
  v <- akash_log_cdf_surv(x, sigma)
  ifelse(v$lS < log(0.5), log1p(-exp(v$lS)), v$lF)
}

# log(-log F(x; sigma)), decreasing in x; equals log S when S is tiny.
akash_lmlF <- function(x, sigma) {
  v <- akash_log_cdf_surv(x, sigma)
  out <- ifelse(v$lS < log(0.5), log(-log1p(-exp(v$lS))), log(-v$lF))
  ifelse(v$lS < -700, v$lS, out)
}

akash_log_pdf <- function(x, sigma) {
  3 * log(sigma) - log(sigma^2 + 2) + log1p(x^2) - sigma * x
}

check_par <- function(mu, sigma, tau) {
  if (any(!is.finite(c(mu, sigma, tau))) || any(mu <= 0) || any(sigma <= 0) ||
      any(tau <= 0 | tau >= 1)) {
    stop("Use mu > 0, sigma > 0 and 0 < tau < 1.")
  }
}

alpha_qeakash <- function(mu, sigma, tau) {
  check_par(mu, sigma, tau)
  log(tau) / akash_log_cdf(mu, sigma)
}

# Density; at x = 0 it returns the one-sided limit (Inf, sigma^3/(sigma^2+2)
# or 0 for alpha < 1, = 1 or > 1).
dQEAkash <- function(x, mu, sigma, tau, log = FALSE) {
  a <- alpha_qeakash(mu, sigma, tau)
  n <- max(length(x), length(a), length(sigma))
  x <- rep_len(x, n); a <- rep_len(a, n); sigma <- rep_len(sigma, n)
  out <- rep(-Inf, n)
  k <- which(x > 0 & is.finite(x))
  out[k] <- log(a[k]) + akash_log_pdf(x[k], sigma[k]) +
    (a[k] - 1) * akash_log_cdf(x[k], sigma[k])
  k0 <- which(x == 0)
  out[k0] <- ifelse(a[k0] < 1, Inf, ifelse(a[k0] > 1, -Inf,
                    3 * log(sigma[k0]) - log(sigma[k0]^2 + 2)))
  if (log) out else exp(out)
}

# CDF G(x) = F(x; sigma)^alpha, or survival 1 - G(x) without cancellation.
pQEAkash <- function(q, mu, sigma, tau, lower.tail = TRUE, log.p = FALSE) {
  a <- alpha_qeakash(mu, sigma, tau)
  lG <- ifelse(q > 0, a * akash_log_cdf(pmax(q, 0), sigma), -Inf)
  out <- if (lower.tail) lG else log(-expm1(lG))
  if (log.p) out else exp(out)
}

hQEAkash <- function(x, mu, sigma, tau) {
  exp(dQEAkash(x, mu, sigma, tau, log = TRUE) -
        pQEAkash(x, mu, sigma, tau, lower.tail = FALSE, log.p = TRUE))
}

# Solve log(-log F(x; sigma)) = target for x by vectorized bisection on
# log(x); start is an initial upper end on the log scale.
akash_solve <- function(target, sigma, start) {
  lo <- rep(-200, length(target))
  hi <- start
  repeat {
    j <- akash_lmlF(exp(lo), sigma) <= target & lo > -740
    if (!any(j)) break
    lo[j] <- pmax(lo[j] * 1.5, -740)
  }
  repeat {
    j <- akash_lmlF(exp(hi), sigma) > target
    if (!any(j)) break
    hi[j] <- hi[j] + 2
  }
  for (it in 1:80) {
    md <- (lo + hi) / 2
    j <- akash_lmlF(exp(md), sigma) > target
    lo[j] <- md[j]
    hi[!j] <- md[!j]
  }
  exp((lo + hi) / 2)
}

# Quantile function Q(p) = F^{-1}(p^(1/alpha); sigma), with no closed form:
# log(-log F(Q(p))) = log(-log p) - log(alpha).
qQEAkash <- function(p, mu, sigma, tau) {
  check_par(mu, sigma, tau)
  n <- max(length(p), length(mu), length(sigma))
  p <- rep_len(p, n); mu <- rep_len(mu, n); s <- rep_len(sigma, n)
  out <- rep(NA_real_, n)
  out[p <= 0] <- 0
  out[p >= 1] <- Inf
  k <- which(p > 0 & p < 1)
  if (!length(k)) return(out)
  target <- log(-log(p[k])) - (log(-log(tau)) - akash_lmlF(mu[k], s[k]))
  out[k] <- akash_solve(target, s[k], log(mu[k]) + 1)
  out
}

rQEAkash <- function(n, mu, sigma, tau) qQEAkash(runif(n), mu, sigma, tau)

# Bowley skewness (quartiles) and Moors kurtosis (octiles).
shape_measures <- function(mu, sigma, tau) {
  q <- qQEAkash((1:7) / 8, mu, sigma, tau)
  c(bowley = (q[6] + q[2] - 2 * q[4]) / (q[6] - q[2]),
    moors = ((q[7] - q[5]) + (q[3] - q[1])) / (q[6] - q[2]))
}

# Raw moment E[X^r] by adaptive quadrature, split at the median. When
# alpha < 1 the substitution v = x^alpha removes the singularity at zero.
moment_qeakash <- function(r, mu, sigma, tau) {
  a <- alpha_qeakash(mu, sigma, tau)
  f <- function(x) x^r * dQEAkash(x, mu, sigma, tau)
  cut <- qQEAkash(0.5, mu, sigma, tau)
  lower <- if (a >= 1) {
    integrate(f, 0, cut, rel.tol = 1e-9, subdivisions = 2000L)$value
  } else {
    integrate(function(v) {
      val <- f(v^(1 / a)) * (1 / a) * v^(1 / a - 1)
      val[!is.finite(val)] <- 0
      val
    }, 0, cut^a, rel.tol = 1e-9, subdivisions = 2000L)$value
  }
  lower + integrate(f, cut, Inf, rel.tol = 1e-9, subdivisions = 2000L)$value
}

# Second derivatives of log g(y; mu, sigma), summed over the sample. With
# L = log F(mu), M = log F(y), a(x) = f/F and b(x) = (dF/dsigma)/F, the score
# is dl/dmu = alpha_mu (1/alpha + M) and dl/dsigma = alpha_sigma (1/alpha + M)
# + 3/sigma - 2 sigma/(sigma^2 + 2) - y + (alpha - 1) b(y).
hessian_qeakash <- function(y, mu, sigma, tau) {
  K <- function(x) (sigma^2 + 2) * (1 + x^2) + 2 * sigma * x + 4
  Ks <- function(x) 2 * sigma * (1 + x^2) + 2 * x
  log_Fs <- function(x) {
    2 * log(sigma) + log(x) - sigma * x + log(K(x)) - 2 * log(sigma^2 + 2)
  }
  L <- akash_log_cdf(mu, sigma)
  M <- akash_log_cdf(y, sigma)
  A <- log(tau) / L
  a_mu <- exp(akash_log_pdf(mu, sigma) - L)
  b_mu <- exp(log_Fs(mu) - L)
  b_y <- exp(log_Fs(y) - M)
  ds <- 3 / sigma - 2 * sigma / (sigma^2 + 2)
  L_mm <- a_mu * (2 * mu / (1 + mu^2) - sigma) - a_mu^2
  L_ms <- a_mu * (ds - mu) - a_mu * b_mu
  L_ss <- b_mu * (2 / sigma - mu + Ks(mu) / K(mu) - 4 * sigma / (sigma^2 + 2)) - b_mu^2
  M_ss <- b_y * (2 / sigma - y + Ks(y) / K(y) - 4 * sigma / (sigma^2 + 2)) - b_y^2
  A_m <- -A * a_mu / L
  A_s <- -A * b_mu / L
  A_mm <- A * (2 * a_mu^2 / L^2 - L_mm / L)
  A_ms <- A * (2 * a_mu * b_mu / L^2 - L_ms / L)
  A_ss <- A * (2 * b_mu^2 / L^2 - L_ss / L)
  V <- 1 / A + M
  U <- c(sum(A_m * V), sum(ds - y + A_s * V + (A - 1) * b_y))
  H <- matrix(c(sum(A_mm * V - A_m^2 / A^2),
                sum(A_ms * V - A_m * A_s / A^2 + A_m * b_y),
                sum(A_ms * V - A_m * A_s / A^2 + A_m * b_y),
                sum(-3 / sigma^2 + 2 * (sigma^2 - 2) / (sigma^2 + 2)^2 + A_ss * V -
                      A_s^2 / A^2 + 2 * A_s * b_y + (A - 1) * M_ss)), 2, 2)
  list(U = U, H = H)
}

loglik_qeakash <- function(y, mu, sigma, tau) {
  sum(dQEAkash(y, mu, sigma, tau, log = TRUE))
}

# Profile log-likelihood in sigma. For fixed sigma, alpha(sigma) = -n/M(sigma)
# with M(sigma) = sum log F(y; sigma); it does not depend on tau.
profile_loglik <- function(sigma, y) {
  n <- length(y)
  vapply(sigma, function(s) {
    M <- sum(akash_log_cdf(y, s))
    if (!is.finite(M) || M >= 0) return(NA_real_)
    n * log(-n / M) + 3 * n * log(s) - n * log(s^2 + 2) + sum(log1p(y^2)) -
      s * sum(y) - n - M
  }, numeric(1))
}

# Maximum likelihood by the profile log-likelihood. sigma maximizes l_p, and
# mu(sigma) = F^{-1}(tau^(1/alpha(sigma)); sigma). Standard errors come from
# the inverse observed information J = -H at the estimate.
fit_qeakash <- function(y, tau) {
  f <- function(ls) {
    v <- profile_loglik(exp(ls), y)
    if (is.finite(v)) -v else 1e100
  }
  grid <- log(1 / mean(y)) + seq(-8, 8, by = 0.5)
  k <- which.min(vapply(grid, f, numeric(1)))
  opt <- optimize(f, grid[k] + c(-0.5, 0.5), tol = 1e-12)
  sigma <- exp(opt$minimum)
  alpha <- -length(y) / sum(akash_log_cdf(y, sigma))
  # mu solves F(mu; sigma) = tau^(1/alpha), that is,
  # log(-log F(mu; sigma)) = log(-log tau) - log(alpha).
  mu <- akash_solve(log(-log(tau)) - log(alpha), sigma, log(max(y)) + 1)
  d <- hessian_qeakash(y, mu, sigma, tau)
  J <- -d$H
  V <- tryCatch(if (all(eigen(J, symmetric = TRUE, only.values = TRUE)$values > 0))
    solve(J) else NULL, error = function(e) NULL)
  list(mu = mu, sigma = sigma, alpha = alpha, loglik = -opt$objective,
       score = d$U, vcov = V)
}

# ============================================================================
# USER INTERFACE
# ============================================================================

presets <- c(
  "Custom" = "custom",
  "Bathtub hazard, alpha < 1 (mu = 0.942, sigma = 1)" = "a",
  "Akash, alpha = 1 (mu = 1.983, sigma = 1)" = "b",
  "Increasing hazard, alpha > 1 (mu = 3.035, sigma = 1)" = "c",
  "Maximum then minimum (mu = 0.169, sigma = 5)" = "d"
)
preset_values <- list(a = c(0.942, 1), b = c(1.983, 1), c = c(3.035, 1),
                      d = c(0.169, 5))

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  titlePanel("QEAkash - Quantile Exponentiated Akash distribution"),
  h5("Distribution, skewness, kurtosis and likelihood"),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.tab_ativa == 'Distribution'",
        h4("Distribution"),
        sliderInput("tau_d", "τ (fixed probability level)", min = 0.01,
                    max = 0.99, value = 0.50, step = 0.01),
        numericInput("mu_d", "μτ (quantile of order τ)", value = 1, min = 0.001,
                     step = 0.05),
        numericInput("sigma_d", "σ (Akash parameter)", value = 1, min = 0.001,
                     step = 0.1),
        selectInput("preset", "Hazard shapes with τ = 0.5", choices = presets),
        checkboxInput("x_auto", "Plot up to the 0.995 quantile", TRUE),
        conditionalPanel(
          condition = "!input.x_auto",
          sliderInput("x_max", "Upper limit of x", min = 0.1, max = 30,
                      value = 5, step = 0.1)
        )
      ),
      conditionalPanel(
        condition = "input.tab_ativa == 'Skewness and Kurtosis'",
        h4("Skewness and kurtosis"),
        radioButtons("curves", "Curves of", inline = TRUE,
                     choices = c("τ" = "tau", "μτ" = "mu")),
        conditionalPanel(
          condition = "input.curves == 'tau'",
          textInput("tau_values", "Values of τ",
                    "0.10, 0.25, 0.50, 0.75, 0.90"),
          numericInput("mu_fixed", "Fixed μτ", value = 0.5, min = 0.001,
                       step = 0.1)
        ),
        conditionalPanel(
          condition = "input.curves == 'mu'",
          textInput("mu_values", "Values of μτ", "0.5, 1.0, 1.5, 2.0, 2.5"),
          sliderInput("tau_fixed", "Fixed τ", min = 0.01, max = 0.99,
                      value = 0.50, step = 0.01)
        ),
        sliderInput("sigma_range", "Range of σ", min = 0.1, max = 20,
                    value = c(2, 10), step = 0.1),
        helpText("For τ = 0.9 and small σ, Moors kurtosis grows quickly",
                 "(about 19.5 at σ = 1 with μτ = 0.5).")
      ),
      conditionalPanel(
        condition = "input.tab_ativa == 'Log-likelihood'",
        h4("Simulation and estimation"),
        sliderInput("tau_l", "τ (fixed)", min = 0.01, max = 0.99, value = 0.50,
                    step = 0.01),
        numericInput("mu_l", "True μτ", 1, min = 0.01, max = 50, step = 0.1),
        numericInput("sigma_l", "True σ", 2, min = 0.01, max = 50, step = 0.1),
        numericInput("n_l", "Sample size", 200, min = 10, max = 5000, step = 10),
        numericInput("seed_l", "Seed", 123, min = 0, max = 2147483646, step = 1),
        actionButton("simulate", "Generate sample and estimate",
                     class = "btn-primary"),
        hr(),
        sliderInput("width_l", "Surface range around the parameters (factor)",
                    1.2, 4, 2, step = 0.1),
        sliderInput("grid_l", "Surface resolution", 25, 80, 45, step = 5),
        helpText("The sample changes only when the button is clicked."),
        downloadButton("sample_csv", "Download sample")
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "tab_ativa", selected = "Distribution",
        tabPanel(
          "Distribution",
          br(),
          tableOutput("dist_summary"),
          fluidRow(
            column(6, plotOutput("dist_dens")),
            column(6, plotOutput("dist_cdf")),
            column(6, plotOutput("dist_surv")),
            column(6, plotOutput("dist_haz"))
          ),
          helpText("α = log(τ)/log F(μτ; σ), where F is the Akash distribution",
                   "function. The dashed lines mark μτ, the point (μτ, τ) on",
                   "the CDF and the hazard asymptote h = σ.")
        ),
        tabPanel(
          "Skewness and Kurtosis",
          plotOutput("bowley"),
          br(),
          plotOutput("moors")
        ),
        tabPanel(
          "Log-likelihood",
          br(),
          helpText("τ is fixed. μτ is the quantile of order τ. The estimates",
                   "maximize the profile log-likelihood in σ, which does not",
                   "depend on τ, with μτ(σ) = F^{-1}(τ^{1/α(σ)}; σ) and",
                   "α(σ) = -n / Σ log F(yᵢ; σ). Standard errors come from the",
                   "observed information; the 95% intervals are built on the",
                   "log scale."),
          verbatimTextOutput("ll_status"),
          tableOutput("ll_estimates"),
          plotOutput("ll_profile", height = "320px"),
          plotlyOutput("ll_surface", height = "530px"),
          plotOutput("ll_contour", height = "400px"),
          sliderInput("hist_bins", "Number of histogram bins", 10, 80, 30,
                      step = 1),
          plotOutput("ll_sample", height = "380px")
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {

  observeEvent(input$preset, {
    if (input$preset != "custom") {
      v <- preset_values[[input$preset]]
      updateSliderInput(session, "tau_d", value = 0.5)
      updateNumericInput(session, "mu_d", value = v[1])
      updateNumericInput(session, "sigma_d", value = v[2])
    }
  })

  # ---------------------------------------------------------------- distribution
  par_d <- reactive({
    validate(need(is.finite(input$mu_d) && input$mu_d > 0, "Use μτ > 0."),
             need(is.finite(input$sigma_d) && input$sigma_d > 0, "Use σ > 0."))
    list(mu = input$mu_d, sigma = input$sigma_d, tau = input$tau_d)
  })

  dd <- reactive({
    p <- par_d()
    a <- tryCatch(alpha_qeakash(p$mu, p$sigma, p$tau), error = function(e) NA)
    validate(need(is.finite(a) && a > 0,
                  "Parameters outside numerical precision; change μτ or σ."))
    x_max <- if (isTRUE(input$x_auto)) qQEAkash(0.995, p$mu, p$sigma, p$tau) else
      input$x_max
    x_min <- min(qQEAkash(1e-4, p$mu, p$sigma, p$tau), x_max / 1000)
    x <- sort(unique(c(exp(seq(log(x_min), log(x_max), length.out = 200)),
                       seq(x_min, x_max, length.out = 300))))
    df <- data.frame(
      x = x,
      dens = dQEAkash(x, p$mu, p$sigma, p$tau),
      cdf = pQEAkash(x, p$mu, p$sigma, p$tau),
      surv = pQEAkash(x, p$mu, p$sigma, p$tau, lower.tail = FALSE),
      haz = hQEAkash(x, p$mu, p$sigma, p$tau)
    )
    list(df = df[is.finite(df$dens) & is.finite(df$haz), ], alpha = a, p = p,
         x_ref = qQEAkash(0.02, p$mu, p$sigma, p$tau))
  })

  # Vertical limit. When alpha < 1 the density and the hazard diverge at the
  # origin, so the limit is set by the curve beyond the 0.02 quantile.
  ylim_div <- function(d, col) {
    v <- d$df[[col]]
    if (d$alpha < 1) c(0, 1.4 * max(v[d$df$x >= d$x_ref])) else c(0, max(v))
  }

  output$dist_summary <- renderTable({
    d <- dd(); p <- d$p; a <- d$alpha
    regime <- if (abs(a - 1) < 1e-8) {
      "α = 1: Akash distribution"
    } else if (a < 1) {
      "α < 1: density and hazard diverge at the origin"
    } else {
      "α > 1: density and hazard start at zero"
    }
    m <- tryCatch(sapply(1:2, moment_qeakash, mu = p$mu, sigma = p$sigma,
                         tau = p$tau), error = function(e) c(NA, NA))
    q <- qQEAkash(c(0.25, 0.5, 0.75), p$mu, p$sigma, p$tau)
    sm <- shape_measures(p$mu, p$sigma, p$tau)
    data.frame(
      Quantity = c("α", "Regime", "G(μτ)", "Mean", "Variance",
                   "Quartiles (Q1, Q2, Q3)", "Bowley skewness",
                   "Moors kurtosis"),
      Value = c(sprintf("%.6g", a), regime,
                sprintf("%.4f", pQEAkash(p$mu, p$mu, p$sigma, p$tau)),
                sprintf("%.6g", m[1]), sprintf("%.6g", m[2] - m[1]^2),
                paste(sprintf("%.4g", q), collapse = ", "),
                sprintf("%.4f", sm[["bowley"]]), sprintf("%.4f", sm[["moors"]])),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)

  output$dist_dens <- renderPlot({
    d <- dd(); df <- d$df
    ggplot(df, aes(x, dens)) +
      geom_line(colour = "#2C3E50", linewidth = 1) +
      geom_vline(xintercept = d$p$mu, linetype = "dashed", colour = "grey50") +
      coord_cartesian(ylim = ylim_div(d, "dens")) +
      labs(title = "Probability density function", x = "x", y = "g(x)") +
      theme_minimal(base_size = 13)
  })

  output$dist_cdf <- renderPlot({
    d <- dd(); df <- d$df; p <- d$p
    ggplot(df, aes(x, cdf)) +
      geom_line(colour = "#18BC9C", linewidth = 1) +
      annotate("segment", x = p$mu, xend = p$mu, y = 0, yend = p$tau,
               linetype = "dashed", colour = "grey50") +
      annotate("segment", x = 0, xend = p$mu, y = p$tau, yend = p$tau,
               linetype = "dashed", colour = "grey50") +
      annotate("point", x = p$mu, y = p$tau, colour = "#E74C3C", size = 3) +
      labs(title = "Cumulative distribution function", x = "x", y = "G(x)") +
      theme_minimal(base_size = 13)
  })

  output$dist_surv <- renderPlot({
    d <- dd()
    ggplot(d$df, aes(x, surv)) +
      geom_line(colour = "#3498DB", linewidth = 1) +
      geom_vline(xintercept = d$p$mu, linetype = "dashed", colour = "grey50") +
      labs(title = "Survival function", x = "x", y = "S(x)") +
      theme_minimal(base_size = 13)
  })

  output$dist_haz <- renderPlot({
    d <- dd(); df <- d$df
    lim <- ylim_div(d, "haz")
    lim[2] <- max(lim[2], d$p$sigma * 1.05)
    ggplot(df, aes(x, haz)) +
      geom_line(colour = "#E74C3C", linewidth = 1) +
      geom_hline(yintercept = d$p$sigma, linetype = "dashed", colour = "grey50") +
      coord_cartesian(ylim = lim) +
      labs(title = "Hazard function", x = "x", y = "h(x)") +
      theme_minimal(base_size = 13)
  })

  # ---------------------------------------------------- skewness and kurtosis
  parse_values <- function(txt) {
    v <- suppressWarnings(as.numeric(strsplit(txt, "[,; ]+")[[1]]))
    unique(v[is.finite(v)])
  }

  shape_data <- reactive({
    s_seq <- seq(input$sigma_range[1], input$sigma_range[2], length.out = 50)
    if (input$curves == "tau") {
      vals <- parse_values(input$tau_values)
      validate(need(length(vals) > 0 && all(vals > 0 & vals < 1),
                    "Enter values of τ in (0, 1), separated by commas."),
               need(is.finite(input$mu_fixed) && input$mu_fixed > 0,
                    "Use a positive μτ."))
      grid <- expand.grid(sigma = s_seq, value = vals)
      grid$mu <- input$mu_fixed; grid$tau <- grid$value
      legend <- "τ"
    } else {
      vals <- parse_values(input$mu_values)
      validate(need(length(vals) > 0 && all(vals > 0),
                    "Enter positive values of μτ, separated by commas."))
      grid <- expand.grid(sigma = s_seq, value = vals)
      grid$mu <- grid$value; grid$tau <- input$tau_fixed
      legend <- "μτ"
    }
    sm <- t(mapply(function(m, s, t) {
      tryCatch(shape_measures(m, s, t), error = function(e) c(NA, NA))
    }, grid$mu, grid$sigma, grid$tau))
    grid$bowley <- sm[, 1]; grid$moors <- sm[, 2]
    grid$curve <- factor(format(grid$value), levels = format(sort(vals)))
    list(df = grid, legend = legend)
  })

  plot_shape <- function(col, ylab) {
    d <- shape_data()
    df <- d$df[is.finite(d$df[[col]]), ]
    req(nrow(df) > 0)
    ggplot(df, aes(sigma, .data[[col]], colour = curve)) +
      geom_line(linewidth = 0.9) +
      scale_colour_viridis_d(option = "mako", begin = 0.1, end = 0.85,
                             name = d$legend) +
      labs(x = "σ", y = ylab) +
      theme_bw(base_size = 13) +
      theme(panel.grid = element_blank())
  }
  output$bowley <- renderPlot(plot_shape("bowley", "Bowley skewness"))
  output$moors <- renderPlot(plot_shape("moors", "Moors kurtosis"))

  # ------------------------------------------------------------ log-likelihood
  experiment <- eventReactive(input$simulate, {
    validate(
      need(is.finite(input$mu_l) && input$mu_l > 0, "Use μτ > 0."),
      need(is.finite(input$sigma_l) && input$sigma_l > 0, "Use σ > 0."),
      need(is.finite(input$n_l) && input$n_l >= 10 && input$n_l <= 5000 &&
             input$n_l == floor(input$n_l), "Use an integer n between 10 and 5000."),
      need(is.finite(input$seed_l) && input$seed_l >= 0 &&
             input$seed_l == floor(input$seed_l), "Use a valid integer seed.")
    )
    withProgress(message = "Simulating and maximizing the log-likelihood",
                 value = 0.2, {
      tryCatch({
        set.seed(input$seed_l)
        y <- rQEAkash(input$n_l, input$mu_l, input$sigma_l, input$tau_l)
        incProgress(0.4)
        fit <- fit_qeakash(y, input$tau_l)
        list(y = y, fit = fit, mu = input$mu_l, sigma = input$sigma_l,
             tau = input$tau_l, seed = input$seed_l)
      }, error = function(e) validate(need(FALSE, conditionMessage(e))))
    })
  }, ignoreNULL = FALSE)

  llgrid <- reactive({
    e <- experiment(); w <- input$width_l
    mus <- exp(seq(log(min(e$mu, e$fit$mu) / w), log(max(e$mu, e$fit$mu) * w),
                   length.out = input$grid_l))
    sigmas <- exp(seq(log(min(e$sigma, e$fit$sigma) / w),
                      log(max(e$sigma, e$fit$sigma) * w),
                      length.out = input$grid_l))
    # plotly expects rows in sigma and columns in mu
    z <- outer(sigmas, mus, Vectorize(function(s, m) {
      tryCatch(loglik_qeakash(e$y, m, s, e$tau), error = function(err) NA_real_)
    }))
    z[!is.finite(z)] <- NA_real_
    list(mu = mus, sigma = sigmas, z = z)
  })

  output$ll_status <- renderText({
    e <- experiment(); g <- llgrid(); f <- e$fit
    text <- sprintf(paste0("Sample: n = %d | seed = %d | τ = %.2f\n",
                           "Maximum log-likelihood = %.6f | α̂ = %.6g | ",
                           "score at the estimate = (%.2e, %.2e)"),
                    length(e$y), e$seed, e$tau, f$loglik, f$alpha,
                    f$score[1], f$score[2])
    if (is.null(f$vcov)) {
      text <- paste(text, "The observed information is not positive definite.",
                    sep = "\n")
    }
    if (any(g$z > f$loglik + 1e-6, na.rm = TRUE)) {
      text <- paste(text, "Warning: the grid has a value above the fitted maximum.",
                    sep = "\n")
    }
    text
  })

  output$ll_estimates <- renderTable({
    e <- experiment(); f <- e$fit
    est <- c(f$mu, f$sigma)
    se <- if (is.null(f$vcov)) c(NA, NA) else sqrt(diag(f$vcov))
    z <- qnorm(0.975)
    data.frame(Parameter = c("μτ", "σ"), True = c(e$mu, e$sigma),
               Estimate = est, SE = se,
               `95% CI (lower)` = est * exp(-z * se / est),
               `95% CI (upper)` = est * exp(z * se / est),
               check.names = FALSE)
  }, digits = 5)

  output$ll_profile <- renderPlot({
    e <- experiment(); f <- e$fit
    s <- exp(seq(log(f$sigma / input$width_l^1.5),
                 log(f$sigma * input$width_l^1.5), length.out = 200))
    df <- data.frame(sigma = s, lp = profile_loglik(s, e$y))
    ggplot(df[is.finite(df$lp), ], aes(sigma, lp)) +
      geom_line(linewidth = 1, colour = "#2C3E50") +
      geom_vline(xintercept = f$sigma, colour = "red") +
      geom_vline(xintercept = e$sigma, colour = "dodgerblue", linetype = "dashed") +
      labs(x = "σ", y = expression(l[p](sigma)),
           title = "Profile log-likelihood (red: estimate, blue: true σ)") +
      theme_minimal(base_size = 13)
  })

  output$ll_surface <- renderPlotly({
    e <- experiment(); g <- llgrid(); f <- e$fit
    true_ll <- loglik_qeakash(e$y, e$mu, e$sigma, e$tau)
    fig <- plot_ly(x = g$mu, y = g$sigma, z = g$z, type = "surface",
                   colorscale = "Viridis", showlegend = FALSE,
                   colorbar = list(title = "ℓ", x = 1.02, y = 0.55, len = 0.8))
    fig <- add_trace(fig, x = f$mu, y = f$sigma, z = f$loglik,
                     type = "scatter3d", mode = "markers", inherit = FALSE,
                     name = "Maximum", showlegend = TRUE,
                     marker = list(color = "red", size = 6))
    fig <- add_trace(fig, x = e$mu, y = e$sigma, z = true_ll, type = "scatter3d",
                     mode = "markers", inherit = FALSE, name = "True parameters",
                     showlegend = TRUE,
                     marker = list(color = "dodgerblue", size = 5,
                                   symbol = "diamond"))
    layout(fig,
           scene = list(xaxis = list(title = "μτ"), yaxis = list(title = "σ"),
                        zaxis = list(title = "ℓ(μτ, σ)")),
           legend = list(orientation = "h", x = 0.5, xanchor = "center",
                         y = -0.12, yanchor = "top"),
           margin = list(l = 10, r = 85, b = 105, t = 25))
  })

  output$ll_contour <- renderPlot({
    e <- experiment(); g <- llgrid()
    df <- expand.grid(sigma = g$sigma, mu = g$mu)
    df$ll <- as.vector(g$z)
    ggplot(df, aes(mu, sigma, z = ll)) +
      geom_contour_filled(bins = 18, na.rm = TRUE) +
      annotate("point", x = e$fit$mu, y = e$fit$sigma, colour = "red", size = 3) +
      annotate("point", x = e$mu, y = e$sigma, colour = "dodgerblue", shape = 18,
               size = 3) +
      labs(x = "μτ", y = "σ", fill = "ℓ", title = "Contour levels") +
      theme_minimal(base_size = 13)
  })

  output$ll_sample <- renderPlot({
    e <- experiment(); f <- e$fit
    x_max <- max(e$y) * 1.05
    h <- hist(e$y, breaks = seq(0, x_max, length.out = input$hist_bins + 1),
              plot = FALSE)
    bars <- data.frame(left = head(h$breaks, -1), right = tail(h$breaks, -1),
                       density = h$density)
    x_min <- max(min(e$y) * 0.1, .Machine$double.xmin)
    x <- sort(unique(c(exp(seq(log(x_min), log(x_max), length.out = 400)),
                       seq(x_min, x_max, length.out = 600))))
    df <- rbind(
      data.frame(x = x, density = dQEAkash(x, f$mu, f$sigma, e$tau),
                 curve = "Fitted"),
      data.frame(x = x, density = dQEAkash(x, e$mu, e$sigma, e$tau),
                 curve = "True"))
    df <- df[is.finite(df$density), ]
    y_max <- 1.25 * max(h$density)
    ggplot() +
      geom_rect(data = bars, aes(xmin = left, xmax = right, ymin = 0,
                                 ymax = density),
                fill = "#B9CFDF", colour = "white", linewidth = 0.3) +
      geom_line(data = df, aes(x, density, colour = curve, linetype = curve),
                linewidth = 1) +
      scale_colour_manual(values = c(Fitted = "#D62728", True = "dodgerblue")) +
      scale_linetype_manual(values = c(Fitted = "solid", True = "dashed")) +
      coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max), expand = FALSE) +
      labs(x = "y", y = "Density", colour = NULL, linetype = NULL,
           title = "Histogram of the simulated sample",
           subtitle = sprintf("Estimates: μτ = %.4f | σ = %.4f | τ = %.2f",
                              f$mu, f$sigma, e$tau)) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })

  output$sample_csv <- downloadHandler(
    filename = function() "QEAkash_sample.csv",
    content = function(file) {
      e <- experiment()
      write.csv(data.frame(y = e$y, mu_tau = e$mu, sigma = e$sigma, tau = e$tau,
                           seed = e$seed), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
