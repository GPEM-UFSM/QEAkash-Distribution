# QEAkash Distribution and QR-QEAkash Regression Model

This repository provides R code and Monte Carlo replication files for the
Quantile Exponentiated Akash (QEAkash) distribution and its quantile
regression model (QR-QEAkash). The distribution is defined for continuous
positive responses. It reparameterizes the exponentiated Akash
distribution by its quantile: for a fixed probability level $\tau$, the
parameter $\mu_\tau$ is the quantile of order $\tau$, and $\sigma$ is the
parameter of the Akash base distribution. With $F(x;\sigma)$ the Akash
distribution function,

$$
G(x;\mu_\tau,\sigma) = F(x;\sigma)^{\alpha},
\qquad
\alpha = \frac{\log(\tau)}{\log F(\mu_\tau;\sigma)},
$$

so that $G(\mu_\tau;\mu_\tau,\sigma) = \tau$. In QR-QEAkash regression,
covariates enter through link functions for the conditional quantile and
for $\sigma$, and right-censored responses are allowed.

## Repository contents

| Directory | Contents |
|:--|:--|
| [Gamlss_Function](Gamlss_Function/) | QEAkash family for GAMLSS and examples of fitting an intercept-only distribution and a conditional quantile regression with submodels for $\mu_\tau$ and $\sigma$. |
| [Moments](Moments/) | R code for computing raw QEAkash moments and comparing adaptive quadrature with the binomial series representation. |
| [Simulations_QEAkash](Simulations_QEAkash/) | R scripts for the QEAkash simulation and post-processing, with 24 archived replication files for scenarios D1 and D2. |
| [Simulations_QR_QEAkash](Simulations_QR_QEAkash/) | R scripts for the QR-QEAkash simulation without censoring and with 15%, 30% and 60% right censoring, with 16 archived replication files for scenarios R1, R2_15, R2_30 and R2_60. |
| [Shiny_App](Shiny_App/) | Interactive plots of the density, cumulative distribution, survival and hazard functions, Bowley skewness and Moors kurtosis, and the log-likelihood of a simulated sample. |
| [Launch Shiny App](https://gpem-ufsm.github.io/QEAkash-Distribution/) | Run the interactive QEAkash application directly in your browser. |

## Getting started

The scripts and application use R. Open the README in the relevant
directory for package requirements and examples.

- **Fit a model:** From [Gamlss_Function](Gamlss_Function/), run
  `Rscript QEAkash_Examples.R`.
- **Compute moments:** From [Moments](Moments/), run
  `Rscript QEAkash_Moments.R` to generate the CSV and LaTeX comparison tables.
- **Fit a censored regression:** Source
  `Simulations_QR_QEAkash/QR_QEAkash_Simulation_Functions/01_FUNCOES_QR_QEAKASH.R`
  and call `qr_qeakash_fit()`.
- **Explore the distribution:** From the repository root, run
  `shiny::runApp("Shiny_App")` in R after installing the app dependencies.

## Monte Carlo replication files

The [QEAkash archive](Simulations_QEAkash/QEAkash_Simulations/) and the
[QR-QEAkash archive](Simulations_QR_QEAkash/QR-QEAkash_Simulations/)
contain one `.rds` file per configuration. Each configuration combines a
scenario, a sample size and a quantile level and contains 10,000
replications, with the estimates of the accepted ones. From the repository
root, use `readRDS()` to inspect an archived configuration:

```r
qeakash <- readRDS("Simulations_QEAkash/QEAkash_Simulations/D1_tau_05_n_050.rds")
qr_qeakash <- readRDS("Simulations_QR_QEAkash/QR-QEAkash_Simulations/R1_tau_05_n_050.rds")
```

Both examples select the first scenario at $\tau=0.50$ and $n=50$. The
simulation READMEs describe the scenario parameters and how to rerun the
studies and generate tables and figures.
