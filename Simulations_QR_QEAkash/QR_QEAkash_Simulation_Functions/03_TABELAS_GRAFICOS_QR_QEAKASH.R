# Tables and figures for the QR-QEAkash regression simulation.
#
# Usage, from this folder:
#   Rscript 03_TABELAS_GRAFICOS_QR_QEAKASH.R [input_dir] [output_dir]
# input_dir defaults to the archived run, ../QR-QEAkash_Simulations, and
# output_dir to TABELAS_E_GRAFICOS/<name of input_dir> in this folder.
#
# For each scenario, sample size and coefficient, the summary reports the
# mean, bias, mean squared error (MSE) and empirical coverage probability
# (CP) of the 95% Wald intervals, computed over the accepted replications,
# the mean censoring proportion and the number of failed replications.
# Outputs:
#   summary_qr_qeakash.csv       one row per scenario, n and coefficient
#   table_qr_qeakash_<R>.tex     LaTeX table body for scenario <R>
#   boxplot_<R>.pdf              estimates by sample size (requires ggplot2)
#   density_<R>.pdf              standardized estimates against N(0, 1)

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
  file.path(script_dir, "..", "QR-QEAkash_Simulations")
input_dir <- normalizePath(input_dir, mustWork = TRUE)
output_dir <- if (length(args) >= 2) args[2] else
  file.path(script_dir, "TABELAS_E_GRAFICOS", basename(input_dir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "^R[0-9].*_n_[0-9]+\\.rds$",
                    full.names = TRUE)
if (!length(files)) stop("No configuration files found in ", input_dir)
runs <- lapply(files, readRDS)
PAR <- c("beta11", "beta12", "beta13", "beta21", "beta22", "beta23")

summary_table <- do.call(rbind, lapply(runs, function(z) {
  e <- z$estimates
  do.call(rbind, lapply(PAR, function(p) {
    truth <- z$truth[[p]]
    data.frame(z$configuration, parameter = p, true = truth,
               mean = mean(e[[p]]), bias = mean(e[[p]]) - truth,
               mse = mean((e[[p]] - truth)^2),
               cp = mean(e[[paste0("covered_", p)]]),
               mean_censoring = mean(e$censoring), accepted = nrow(e),
               failures = z$failures)
  }))
}))
summary_table <- summary_table[order(summary_table$scenario,
                                     match(summary_table$parameter, PAR),
                                     summary_table$n), ]
utils::write.csv(summary_table, file.path(output_dir, "summary_qr_qeakash.csv"),
                 row.names = FALSE)
options(width = 160)
print(format(summary_table, digits = 4), row.names = FALSE)

# LaTeX rows: four decimals, with scientific notation for an MSE below 5e-4
# and for a bias that would round to zero.
fmt <- function(x, below = 0) {
  if (x == 0 || abs(x) >= below) return(formatC(x, format = "f", digits = 4))
  e <- floor(log10(abs(x)))
  sprintf("%s\\times10^{%d}", formatC(x / 10^e, format = "f", digits = 1), e)
}
for (scenario in unique(summary_table$scenario)) {
  d <- summary_table[summary_table$scenario == scenario, ]
  lines <- c("\\begin{tabular}{cl rrrrr}", "\\toprule",
             "Parameter & $n$ & Mean & Bias & MSE & CP & Failures \\\\",
             "\\midrule")
  for (p in PAR) {
    dp <- d[d$parameter == p, ]
    if (p != PAR[1]) lines <- c(lines, "\\midrule")
    for (i in seq_len(nrow(dp))) {
      r <- dp[i, ]
      lines <- c(lines, sprintf(
        "%s & $%d$ & $%s$ & $%s$ & $%s$ & $%s$ & $%d$ \\\\",
        if (i == 1) sprintf("\\multirow{%d}{*}{$\\beta_{%s}$}", nrow(dp),
                            substr(p, 5, 6)) else "",
        r$n, fmt(r$mean), fmt(r$bias, 5e-5), fmt(r$mse, 5e-4), fmt(r$cp),
        r$failures))
    }
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file.path(output_dir,
                              sprintf("table_qr_qeakash_%s.tex", scenario)))
}

# Figures: boxplots of the estimates (without outliers) with the true value
# as a dashed line and the mean as a white diamond, and kernel densities of
# the estimates standardized within each configuration against N(0, 1).
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  labels <- setNames(sprintf("beta[%s]", substr(PAR, 5, 6)), PAR)
  for (scenario in unique(summary_table$scenario)) {
    sel <- Filter(function(z) z$configuration$scenario == scenario, runs)
    long <- do.call(rbind, lapply(sel, function(z) {
      e <- z$estimates
      n <- factor(z$configuration$n, levels = sort(unique(summary_table$n)))
      do.call(rbind, lapply(PAR, function(p) {
        data.frame(parameter = factor(labels[[p]], levels = labels), n = n,
                   estimate = e[[p]],
                   z = (e[[p]] - mean(e[[p]])) / sd(e[[p]]),
                   truth = z$truth[[p]])
      }))
    }))
    truths <- unique(long[c("parameter", "truth")])
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
      facet_wrap(~ parameter, scales = "free_x", labeller = label_parsed,
                 nrow = 2) +
      labs(x = "Estimate", y = "n") + theme_bw()
    ggsave(file.path(output_dir, paste0("boxplot_", scenario, ".pdf")), p,
           width = 9, height = 5.5)
    p <- ggplot(long, aes(x = z, colour = n)) +
      geom_density() +
      stat_function(fun = stats::dnorm, colour = "black", linetype = "dashed") +
      facet_wrap(~ parameter, labeller = label_parsed, nrow = 2) +
      labs(x = "Standardized estimate", y = "Density", colour = "n") +
      theme_bw()
    ggsave(file.path(output_dir, paste0("density_", scenario, ".pdf")), p,
           width = 9, height = 5.5)
  }
} else {
  message("ggplot2 is not installed; figures were skipped.")
}
cat("Tables and figures written to", output_dir, "\n")
