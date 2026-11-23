# QEAkash distribution: Monte Carlo simulation scenarios

The simulation evaluates maximum likelihood estimation of the quantile
parameter `mu_tau` and of the Akash parameter `sigma` in the Quantile
Exponentiated Akash (QEAkash) distribution. The probability level `tau` is
fixed within each configuration.

## Simulation design

| Scenario | `sigma` | `mu_tau` at `tau = 0.1, 0.5, 0.9` | `alpha` at `tau = 0.1, 0.5, 0.9` | Configurations | Valid replications |
|:--|--:|:--|:--|--:|--:|
| D1 | 2.98 | 1.25, 1.79, 2.60 | 31.10, 30.92, 30.69 | 12 | 120,000 |
| D2 | 10.08 | 0.05, 0.14, 0.32 | 2.44, 2.35, 2.23 | 12 | 120,000 |
| **Total** | | | | **24** | **240,000** |

Within each scenario, all four sample sizes (`n = 50, 100, 200, 400`) were
crossed with the three quantile levels. Each combination of scenario,
quantile level and sample size defines one configuration, with 10,000
replications. In each scenario, the three `mu_tau` values are quantiles of
one distribution, so `alpha = log(tau)/log F(mu_tau; sigma)` stays close
across levels.

## Data generation and estimation

- For each replication, an independent sample of size `n` was drawn by
  inverse transform sampling, `X = Q(U)` with `U ~ U(0, 1)`. The quantile
  function solves `log(-log F(x; sigma)) = log(-log U) - log(alpha)` by
  bisection on `log(x)`, with the Akash distribution function computed as
  a mixture of `Exp(sigma)` and `Gamma(3, sigma)`.
- The log-likelihood was maximized in `(log mu_tau, log sigma)` with
  `nlminb()` and the closed-form score. The starting value of `mu_tau` is
  the sample quantile of order `tau`, and that of `sigma` comes from a
  one-dimensional profile search. Up to three Newton steps with the
  observed information follow, and the maximum is accepted when the Newton
  decrement is below `1e-6` and the observed information is positive
  definite. Otherwise the replication counts as a failure.
- Standard errors come from the inverse observed information, and the 95%
  Wald intervals are built for `log(mu_tau)` and `log(sigma)`.
- The seed of replication `k` is
  `60000 + 1000 id + 100 k_tau + n + 10 k`, where `id` is 1 for D1 and 2
  for D2 and `k_tau` is 1, 2 or 3 for `tau = 0.1, 0.5, 0.9`. No replication
  failed in the archived run.

## Reported measures

For each configuration, the study summarizes the mean, bias, mean squared
error (MSE) and empirical coverage probability (CP) of the 95% Wald
intervals, for the estimates of `mu_tau` and `sigma`, and the number of
failed replications.

## Reproduction

The simulation depends on R and the `numDeriv` package, and the figures on
`ggplot2`. Install them with `install.packages(c("numDeriv", "ggplot2"))`.
The scripts were tested with R 4.5.3.

| Script | Purpose |
|:-------|:--------|
| `01_FUNCOES_QEAKASH.R` | Stable Akash distribution function, data generation, log-likelihood, score, maximum likelihood estimation and one replication |
| `02_SIMULACAO_QEAKASH.R` | Scenarios, seeded replications in parallel, checkpoints and run log |
| `03_TABELAS_GRAFICOS_QEAKASH.R` | CSV summary, LaTeX tables, boxplots and standardized density plots |

Run the following commands from
[QEAkash_Simulation_Functions](QEAkash_Simulation_Functions/). The first
runs a small test with 40 replications of D1 at `tau = 0.5` and `n = 50`.
The second generates the tables and figures for that test.

```sh
Rscript -e 'options(QEAKASH.mode="test", QEAKASH.scenarios="D1", QEAKASH.tau=0.5, QEAKASH.n=50); source("02_SIMULACAO_QEAKASH.R")'
Rscript 03_TABELAS_GRAFICOS_QEAKASH.R RESULTADOS_MONTE_CARLO/test_qeakash_R00040
```

Because the seeds depend only on the configuration and the replication
index, these 40 replications reproduce the first 40 of the archived run.
Run the full 24-configuration study and generate its tables and figures
with the following commands. The full run takes a few minutes on a
many-core machine and longer on a laptop.

```sh
Rscript 02_SIMULACAO_QEAKASH.R
Rscript 03_TABELAS_GRAFICOS_QEAKASH.R RESULTADOS_MONTE_CARLO/full_qeakash_R10000
```

Rerun the same simulation command to resume an interrupted run. Set the R
option `QEAKASH.workers` to control the number of parallel processes and
`QEAKASH.output_root` to direct a run to another directory. Without
arguments, `03_TABELAS_GRAFICOS_QEAKASH.R` summarizes the archived run and
writes to `TABELAS_E_GRAFICOS/QEAkash_Simulations/`.

## Archived results

The [archive](QEAkash_Simulations/) contains one `.rds` file per
configuration, named `<scenario>_tau_<tau>_n_<n>.rds`. For example,
`D1_tau_05_n_050.rds` is D1 at `tau = 0.5` and `n = 50`. Each file is a
list with the following elements:

| Element | Contents |
|:--|:--|
| `configuration` | Scenario, `tau` and `n` |
| `truth` | True `mu_tau` and `sigma` |
| `alpha` | True power parameter |
| `replications`, `failures` | Number of replications and of failed fits |
| `seeds` | Seed formula |
| `estimates` | One row per accepted replication: index, seed, estimates of `mu_tau` and `sigma`, and whether each 95% Wald interval covered the true value |

```r
d1 <- readRDS("QEAkash_Simulations/D1_tau_05_n_050.rds")
colMeans(d1$estimates[c("mu_tau", "sigma")])
mean(d1$estimates$covered_mu_tau)
```
