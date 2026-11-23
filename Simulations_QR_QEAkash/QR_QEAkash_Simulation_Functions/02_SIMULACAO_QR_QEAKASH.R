# Monte Carlo simulation of maximum likelihood estimation in the QR-QEAkash
# regression model, with tau = 0.5 and logarithmic links in both submodels.
# Each configuration combines a scenario and a sample size n and is saved
# as one .rds file. Completed files are skipped, so rerunning the same
# command resumes an interrupted run.
#
# Settings are read from R options (defaults in parentheses):
#   QRQEAKASH.mode          "full" (10,000 replications) or "test" (40)
#   QRQEAKASH.scenarios     scenarios to run ("R1", "R2_15", "R2_30", "R2_60")
#   QRQEAKASH.n             sample sizes (50, 100, 200, 400)
#   QRQEAKASH.replications  overrides the number of replications
#   QRQEAKASH.workers       parallel R processes (all cores but one)
#   QRQEAKASH.output_root   output directory (RESULTADOS_MONTE_CARLO here)
#
# Example, from this folder:
#   Rscript -e 'options(QRQEAKASH.mode="test", QRQEAKASH.scenarios="R1"); source("02_SIMULACAO_QR_QEAKASH.R")'

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
source(file.path(script_dir, "01_FUNCOES_QR_QEAKASH.R"))
for (p in c("numDeriv", "survival")) if (!requireNamespace(p, quietly = TRUE)) {
  stop("Install ", p, " with install.packages(\"", p, "\").")
}

# Covariate pools. R1 uses the 113 hospitals of the SENIC study (Kutner et
# al., Applied Linear Statistical Models, 5th ed., Appendix C1): average age
# and infection risk in log(mu), and indicators of the North-Central and
# South regions in log(sigma). The R2 scenarios use the first 312 patients
# of survival::pbc (Mayo Clinic primary biliary cholangitis trial): age and
# log(bilirubin) in log(mu), and log(bilirubin) and log(albumin) in log(sigma).
senic <- utils::read.table(file.path(script_dir, "senic.txt"),
                           comment.char = "#", skip = 16)
names(senic) <- c("hospital", "stay", "age", "infection", "culture", "xray",
                  "beds", "school", "region", "patients", "nurses",
                  "facilities")
stopifnot(nrow(senic) == 113)
pbc <- survival::pbc[1:312, ]

pbc_scenario <- function(id, censoring) {
  list(id = id, censoring = censoring,
       X1 = cbind(1, pbc$age, log(pbc$bili)),
       X2 = cbind(1, log(pbc$bili), log(pbc$albumin)),
       beta1 = c(8.81, -0.01, -0.51), beta2 = c(-4.48, 0.59, -2.41))
}
# The id fixes the random seeds. All scenarios share tau = 0.5.
SCENARIOS <- list(
  R1 = list(id = 1, censoring = 0,
            X1 = cbind(1, senic$age, senic$infection),
            X2 = cbind(1, senic$region == 2, senic$region == 3),
            beta1 = c(1.47, 0.01, 0.05), beta2 = c(-0.01, 0.31, 0.47)),
  R2_15 = pbc_scenario(id = 3, censoring = 0.15),
  R2_30 = pbc_scenario(id = 4, censoring = 0.30),
  R2_60 = pbc_scenario(id = 2, censoring = 0.60)
)
NS <- c(50, 100, 200, 400)
TAU <- 0.5
PAR <- c("beta11", "beta12", "beta13", "beta21", "beta22", "beta23")

opt <- function(name, default) getOption(paste0("QRQEAKASH.", name), default)
mode <- match.arg(opt("mode", "full"), c("full", "test"))
scenarios <- opt("scenarios", names(SCENARIOS))
ns <- opt("n", NS)
B <- opt("replications", if (mode == "test") 40L else 10000L)
workers <- opt("workers", max(1L, parallel::detectCores() - 1L))
output_root <- opt("output_root", file.path(script_dir, "RESULTADOS_MONTE_CARLO"))
stopifnot(all(scenarios %in% names(SCENARIOS)), all(ns %in% NS), B >= 1)
out_dir <- file.path(output_root, sprintf("%s_qr_qeakash_R%05d", mode, B))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Output:", out_dir, "\nReplications per configuration:", B,
    "| workers:", workers, "\n")
cl <- NULL
if (workers > 1) {
  cl <- parallel::makeCluster(workers)
  parallel::clusterCall(cl, function(f) source(f),
                        file.path(script_dir, "01_FUNCOES_QR_QEAKASH.R"))
}
# The formal argument starts with a dot so that named arguments passed
# through the dots cannot partially match it.
apply_seeds <- function(.seeds, ...) {
  if (is.null(cl)) lapply(.seeds, qe_replicate_regression, ...) else
    parallel::parLapply(cl, .seeds, qe_replicate_regression, ...)
}

log_rows <- list()
for (scenario in scenarios) {
  sc <- SCENARIOS[[scenario]]
  nu <- qe_calibrate_nu(sc, sc$censoring, tau = TAU)
  cat(sprintf("%s: nu = %s\n", scenario,
              if (is.finite(nu)) format(round(nu, 2)) else "Inf (no censoring)"))
  for (n in ns) {
    destination <- file.path(out_dir, sprintf("%s_tau_05_n_%03d.rds", scenario, n))
    if (file.exists(destination)) {
      cat("Already completed, skipping:", basename(destination), "\n")
      next
    }
    t0 <- Sys.time()
    # Seed of replication k: 90000 + 1000 id + n + 10 k.
    seeds <- 90000 + 1000 * sc$id + n + 10 * seq_len(B)
    out <- do.call(rbind, apply_seeds(seeds, scenario = sc, nu = nu, n = n,
                                      tau = TAU))
    ok <- stats::complete.cases(out)
    estimates <- data.frame(replication = which(ok), seed = seeds[ok],
                            out[ok, 1:6, drop = FALSE],
                            out[ok, 7:12, drop = FALSE] == 1,
                            censoring = out[ok, 13])
    names(estimates)[3:14] <- c(PAR, paste0("covered_", PAR))
    result <- list(
      configuration = data.frame(scenario = scenario, tau = TAU, n = n,
                                 censoring_target = sc$censoring),
      truth = setNames(c(sc$beta1, sc$beta2), PAR),
      censoring_nu = nu,
      replications = B,
      failures = sum(!ok),
      seeds = "90000 + 1000 id + n + 10 k, k = 1, ..., replications",
      estimates = estimates
    )
    saveRDS(result, destination)
    minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("%s: %d failures, mean censoring %.1f%%, %.1f min\n",
                basename(destination), result$failures,
                100 * mean(estimates$censoring), minutes))
    log_rows[[length(log_rows) + 1]] <- data.frame(
      scenario = scenario, n = n, replications = B, failures = result$failures,
      mean_censoring = mean(estimates$censoring), minutes = minutes)
  }
}
if (length(log_rows)) {
  log_file <- file.path(out_dir, "run_log.csv")
  log <- do.call(rbind, log_rows)
  if (file.exists(log_file)) log <- rbind(utils::read.csv(log_file), log)
  utils::write.csv(log, log_file, row.names = FALSE)
}
if (!is.null(cl)) parallel::stopCluster(cl)
cat("Done.\n")
