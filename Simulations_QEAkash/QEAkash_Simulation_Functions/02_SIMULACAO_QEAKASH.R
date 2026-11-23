# Monte Carlo simulation of maximum likelihood estimation in the QEAkash
# distribution. Each configuration crosses a scenario, a quantile level tau
# and a sample size n, and is saved as one .rds file. Completed files are
# skipped, so rerunning the same command resumes an interrupted run.
#
# Settings are read from R options (defaults in parentheses):
#   QEAKASH.mode          "full" (10,000 replications) or "test" (40)
#   QEAKASH.scenarios     scenarios to run ("D1", "D2")
#   QEAKASH.tau           quantile levels (0.1, 0.5, 0.9)
#   QEAKASH.n             sample sizes (50, 100, 200, 400)
#   QEAKASH.replications  overrides the number of replications
#   QEAKASH.workers       parallel R processes (all cores but one)
#   QEAKASH.output_root   output directory (RESULTADOS_MONTE_CARLO here)
#
# Example, from this folder:
#   Rscript -e 'options(QEAKASH.mode="test", QEAKASH.scenarios="D1"); source("02_SIMULACAO_QEAKASH.R")'

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
source(file.path(script_dir, "01_FUNCOES_QEAKASH.R"))
if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("Install numDeriv with install.packages(\"numDeriv\").")
}

# Scenarios: sigma and the true quantile mu_tau at each level tau. The id
# fixes the random seeds.
SCENARIOS <- list(
  D1 = list(id = 1, sigma = 2.98,
            mu = c("0.1" = 1.25, "0.5" = 1.79, "0.9" = 2.60)),
  D2 = list(id = 2, sigma = 10.08,
            mu = c("0.1" = 0.05, "0.5" = 0.14, "0.9" = 0.32))
)
TAUS <- c(0.1, 0.5, 0.9)
NS <- c(50, 100, 200, 400)

opt <- function(name, default) getOption(paste0("QEAKASH.", name), default)
mode <- match.arg(opt("mode", "full"), c("full", "test"))
scenarios <- opt("scenarios", names(SCENARIOS))
taus <- opt("tau", TAUS)
ns <- opt("n", NS)
B <- opt("replications", if (mode == "test") 40L else 10000L)
workers <- opt("workers", max(1L, parallel::detectCores() - 1L))
output_root <- opt("output_root", file.path(script_dir, "RESULTADOS_MONTE_CARLO"))
stopifnot(all(scenarios %in% names(SCENARIOS)), all(taus %in% TAUS),
          all(ns %in% NS), B >= 1)
out_dir <- file.path(output_root, sprintf("%s_qeakash_R%05d", mode, B))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tau_code <- function(tau) sub("^0\\.", "0", format(tau))
config_file <- function(scenario, tau, n) {
  file.path(out_dir, sprintf("%s_tau_%s_n_%03d.rds", scenario, tau_code(tau), n))
}

cat("Output:", out_dir, "\nReplications per configuration:", B,
    "| workers:", workers, "\n")
cl <- NULL
if (workers > 1) {
  cl <- parallel::makeCluster(workers)
  parallel::clusterCall(cl, function(f) source(f),
                        file.path(script_dir, "01_FUNCOES_QEAKASH.R"))
}
# The formal argument starts with a dot so that "s = ..." passed through the
# dots cannot partially match it.
apply_seeds <- function(.seeds, ...) {
  if (is.null(cl)) lapply(.seeds, qe_replicate_distribution, ...) else
    parallel::parLapply(cl, .seeds, qe_replicate_distribution, ...)
}

log_rows <- list()
for (scenario in scenarios) for (tau in taus) for (n in ns) {
  sc <- SCENARIOS[[scenario]]
  mu <- sc$mu[[as.character(tau)]]
  destination <- config_file(scenario, tau, n)
  if (file.exists(destination)) {
    cat("Already completed, skipping:", basename(destination), "\n")
    next
  }
  t0 <- Sys.time()
  # Seed of replication k: 60000 + 1000 id + 100 (index of tau) + n + 10 k.
  seeds <- 60000 + 1000 * sc$id + 100 * match(tau, TAUS) + n + 10 * seq_len(B)
  out <- do.call(rbind, apply_seeds(seeds, mu = mu, s = sc$sigma, tau = tau,
                                    n = n))
  ok <- stats::complete.cases(out)
  result <- list(
    configuration = data.frame(scenario = scenario, tau = tau, n = n),
    truth = c(mu_tau = mu, sigma = sc$sigma),
    alpha = qe_alpha(mu, sc$sigma, tau),
    replications = B,
    failures = sum(!ok),
    seeds = "60000 + 1000 id + 100 k_tau + n + 10 k, k = 1, ..., replications",
    estimates = data.frame(replication = which(ok), seed = seeds[ok],
                           mu_tau = out[ok, 1], sigma = out[ok, 2],
                           covered_mu_tau = out[ok, 3] == 1,
                           covered_sigma = out[ok, 4] == 1)
  )
  saveRDS(result, destination)
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("%s: %d failures, %.1f min\n", basename(destination),
              result$failures, minutes))
  log_rows[[length(log_rows) + 1]] <- data.frame(
    scenario = scenario, tau = tau, n = n, replications = B,
    failures = result$failures, minutes = minutes)
}
if (length(log_rows)) {
  log_file <- file.path(out_dir, "run_log.csv")
  log <- do.call(rbind, log_rows)
  if (file.exists(log_file)) log <- rbind(utils::read.csv(log_file), log)
  utils::write.csv(log, log_file, row.names = FALSE)
}
if (!is.null(cl)) parallel::stopCluster(cl)
cat("Done.\n")
