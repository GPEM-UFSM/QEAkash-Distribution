# Tables and figures for the QEAkash distribution simulation.
#
# Usage, from this folder:
#   Rscript 03_TABELAS_GRAFICOS_QEAKASH.R [input_dir] [output_dir]
# input_dir defaults to the archived run, ../QEAkash_Simulations, and
# output_dir to TABELAS_E_GRAFICOS/<name of input_dir> in this folder.
#
# For each configuration, the summary reports the mean, bias, mean squared
# error (MSE) and empirical coverage probability (CP) of the 95% Wald
# intervals for log(mu_tau) and log(sigma), computed over the accepted
# replications, and the number of failed replications. Outputs:
#   summary_qeakash.csv       one row per configuration
#   table_qeakash_<D>.tex     LaTeX table body for scenario <D>
#   boxplot_<D>_tau_<t>.pdf   estimates by sample size (requires ggplot2)
#   density_<D>_tau_<t>.pdf   standardized estimates against N(0, 1)

script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) {
    return(dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", f[1]),
                                      fixed = TRUE))))
  }
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  normalizePath(getwd())
})()
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[1] else
  file.path(script_dir, "..", "QEAkash_Simulations")
input_dir <- normalizePath(input_dir, mustWork = TRUE)
output_dir <- if (length(args) >= 2) args[2] else
  file.path(script_dir, "TABELAS_E_GRAFICOS", basename(input_dir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "^D[0-9]+_tau_.*\\.rds$",
                    full.names = TRUE)
if (!length(files)) stop("No configuration files found in ", input_dir)
runs <- lapply(files, readRDS)

summary_rows <- lapply(runs, function(z) {
  e <- z$estimates
  measures <- function(est, covered, truth) {
    c(mean = mean(est), bias = mean(est) - truth,
      mse = mean((est - truth)^2), cp = mean(covered))
  }
  m <- measures(e$mu_tau, e$covered_mu_tau, z$truth[["mu_tau"]])
  s <- measures(e$sigma, e$covered_sigma, z$truth[["sigma"]])
  data.frame(z$configuration, mu_tau_true = z$truth[["mu_tau"]],
             sigma_true = z$truth[["sigma"]], alpha_true = z$alpha,
             t(setNames(m, paste0("mu_tau_", names(m)))),
             t(setNames(s, paste0("sigma_", names(s)))),
             accepted = nrow(e), failures = z$failures)
})
summary_table <- do.call(rbind, summary_rows)
summary_table <- summary_table[order(summary_table$scenario, summary_table$tau,
                                     summary_table$n), ]
utils::write.csv(summary_table, file.path(output_dir, "summary_qeakash.csv"),
                 row.names = FALSE)
options(width = 160)
print(format(summary_table, digits = 4), row.names = FALSE)

# LaTeX rows: four decimals, with scientific notation for an MSE below 5e-4
# and for a bias that would round to zero.
fmt_sci <- function(x, below) {
  if (x == 0 || abs(x) >= below) return(formatC(x, format = "f", digits = 4))
  e <- floor(log10(abs(x)))
  sprintf("%s\\times10^{%d}", formatC(x / 10^e, format = "f", digits = 1), e)
}
fmt <- function(x) formatC(x, format = "f", digits = 4)
fmt_bias <- function(x) fmt_sci(x, 5e-5)
fmt_mse <- function(x) fmt_sci(x, 5e-4)
for (scenario in unique(summary_table$scenario)) {
  d <- summary_table[summary_table$scenario == scenario, ]
  lines <- c(
    "\\begin{tabular}{cl*{4}{r}*{4}{r}r}", "\\toprule",
    paste("& & \\multicolumn{4}{c}{$\\widehat{\\mu}_\\tau$} &",
          "\\multicolumn{4}{c}{$\\widehat{\\sigma}$} & \\\\"),
    "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
    paste("$\\tau$ & $n$ & Mean & Bias & MSE & CP & Mean & Bias & MSE & CP &",
          "Failures \\\\"), "\\midrule")
  for (tau in unique(d$tau)) {
    dt <- d[d$tau == tau, ]
    if (tau != d$tau[1]) lines <- c(lines, "\\midrule")
    for (i in seq_len(nrow(dt))) {
      r <- dt[i, ]
      lines <- c(lines, sprintf(
        "%s & $%d$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%d$ \\\\",
        if (i == 1) sprintf("\\multirow{%d}{*}{$%s$}", nrow(dt), format(tau)) else "",
        r$n, fmt(r$mu_tau_mean), fmt_bias(r$mu_tau_bias), fmt_mse(r$mu_tau_mse),
        fmt(r$mu_tau_cp), fmt(r$sigma_mean), fmt_bias(r$sigma_bias),
        fmt_mse(r$sigma_mse), fmt(r$sigma_cp), r$failures))
    }
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file.path(output_dir, sprintf("table_qeakash_%s.tex", scenario)))
}

# Figures: boxplots of the estimates (without outliers) with the true value
# as a dashed line and the mean as a white diamond, and kernel densities of
# the estimates standardized within each configuration against N(0, 1).
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  key <- unique(summary_table[c("scenario", "tau")])
  for (i in seq_len(nrow(key))) {
    sel <- Filter(function(z) z$configuration$scenario == key$scenario[i] &&
                    z$configuration$tau == key$tau[i], runs)
    long <- do.call(rbind, lapply(sel, function(z) {
      e <- z$estimates
      n <- factor(z$configuration$n, levels = sort(unique(summary_table$n)))
      rbind(data.frame(parameter = "mu[tau]", n = n, estimate = e$mu_tau,
                       z = (e$mu_tau - mean(e$mu_tau)) / sd(e$mu_tau),
                       truth = z$truth[["mu_tau"]]),
            data.frame(parameter = "sigma", n = n, estimate = e$sigma,
                       z = (e$sigma - mean(e$sigma)) / sd(e$sigma),
                       truth = z$truth[["sigma"]]))
    }))
    truths <- unique(long[c("parameter", "truth")])
    tag <- sprintf("%s_tau_%s", key$scenario[i], sub("^0\\.", "0", format(key$tau[i])))
    # Box statistics are computed beforehand so that the axes span the
    # whiskers, not the outliers.
    boxes <- do.call(rbind, lapply(split(long, list(long$parameter, long$n),
                                         drop = TRUE), function(g) {
      st <- grDevices::boxplot.stats(g$estimate)$stats
      data.frame(parameter = g$parameter[1], n = g$n[1], xmin = st[1],
                 xlower = st[2], xmiddle = st[3], xupper = st[4], xmax = st[5],
                 mean = mean(g$estimate))
    }))
    p <- ggplot(boxes, aes(y = n)) +
      geom_boxplot(aes(xmin = xmin, xlower = xlower, xmiddle = xmiddle,
                       xupper = xupper, xmax = xmax), stat = "identity",
                   orientation = "y", fill = "grey85", width = 0.6) +
      geom_point(aes(x = mean), shape = 23, fill = "white") +
      geom_vline(data = truths, aes(xintercept = truth), linetype = "dashed",
                 colour = "firebrick") +
      facet_wrap(~ parameter, scales = "free_x", labeller = label_parsed) +
      labs(x = "Estimate", y = "n") + theme_bw()
    ggsave(file.path(output_dir, paste0("boxplot_", tag, ".pdf")), p,
           width = 8, height = 3.5)
    p <- ggplot(long, aes(x = z, colour = n)) +
      geom_density() +
      stat_function(fun = stats::dnorm, colour = "black", linetype = "dashed") +
      facet_wrap(~ parameter, labeller = label_parsed) +
      labs(x = "Standardized estimate", y = "Density", colour = "n") +
      theme_bw()
    ggsave(file.path(output_dir, paste0("density_", tag, ".pdf")), p,
           width = 8, height = 3.5)
  }
} else {
  message("ggplot2 is not installed; figures were skipped.")
}
cat("Tables and figures written to", output_dir, "\n")
