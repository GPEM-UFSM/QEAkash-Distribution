# QEAkash Shiny App

This R Shiny app explores the Quantile Exponentiated Akash (QEAkash)
distribution on the positive real line. Use it to compare the density,
distribution, survival and hazard curves for selected parameters, to
visualize quantile-based shape measures, and to simulate a sample and
inspect its log-likelihood and maximum likelihood estimates.

## Run locally

Install the required R packages once:

```r
install.packages(c("shiny", "bslib", "ggplot2", "plotly"))
```

From the directory containing [app.R](app.R), start the app:

```r
shiny::runApp(".")
```

If your working directory is the repository root, use
`shiny::runApp("Shiny_App")`.

## Parameters and controls

The parameter `mu_tau` is the quantile at probability level `tau`, so
`G(mu_tau) = tau`. For example, `tau = 0.50` makes `mu_tau` the median. The
parameter `sigma` is the parameter of the Akash base distribution, and the
power of the exponentiated Akash distribution is
`alpha = log(tau) / log F(mu_tau; sigma)`, where `F` is the Akash
distribution function.

| Tab | Controls | Output |
|:--|:--|:--|
| **Distribution** | `tau`, `mu_tau`, `sigma`, four hazard-shape presets with `tau = 0.5`, and the plotted range of `x`. | `alpha`, its regime, mean, variance, quartiles, Bowley skewness, Moors kurtosis, and the density, cumulative distribution, survival and hazard functions. |
| **Skewness and Kurtosis** | Curves of `tau` with `mu_tau` fixed, or curves of `mu_tau` with `tau` fixed, and the range of `sigma`. | Bowley skewness and Moors kurtosis as functions of `sigma`. |
| **Log-likelihood** | True `tau`, `mu_tau` and `sigma`, sample size, seed, and the range and resolution of the surface. | Estimates with standard errors and 95% intervals, the profile log-likelihood in `sigma`, the log-likelihood surface and contours, and the histogram of the sample with the fitted and true densities. |

The value of `alpha` governs the behavior at the origin. When
`alpha < 1`, that is, `mu_tau < F^{-1}(tau; sigma)`, the density and the
hazard diverge at zero. When `alpha = 1`, the distribution reduces to the
Akash distribution. When `alpha > 1`, both start at zero. The hazard
function always approaches `sigma` from below as `x` grows.

The **Log-likelihood** tab maximizes the profile log-likelihood
`l_p(sigma)`, obtained by replacing `alpha` with its maximizer
`-n / sum(log F(y_i; sigma))`. It does not depend on `tau`, so changing
`tau` changes the estimated quantile but not the fitted distribution. The
standard errors come from the inverse observed information, computed with
the closed-form second derivatives of the log-likelihood, and the 95%
intervals are built on the log scale.

### Example

In **Distribution**, choose the preset *Akash, alpha = 1* (`tau = 0.50`,
`mu_tau = 1.983`, `sigma = 1`). The cumulative distribution curve passes
through `G(1.983) = 0.50`, and the hazard reaches its minimum at
`sqrt(2) - 1`, about 0.414. Lower `mu_tau` to 0.942 to obtain a bathtub
hazard that diverges at the origin, or choose the last preset to see a
hazard with a local maximum followed by a local minimum.

In **Skewness and Kurtosis**, the default curves of `tau` with
`mu_tau = 0.5` show that both measures approach common values as `sigma`
grows, with Bowley skewness near 0.12 and Moors kurtosis near 1.28.
