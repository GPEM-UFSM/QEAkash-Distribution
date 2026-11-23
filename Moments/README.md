# Raw moments of the QEAkash distribution

[QEAkash_Moments.R](QEAkash_Moments.R) computes raw moments
$\mu_r' = \mathrm{E}[X^r]$ of integer order $r \geq 0$ for the Quantile
Exponentiated Akash (QEAkash) distribution on $(0,\infty)$, by adaptive
quadrature and by a binomial series. It uses base R, plus the optional
`Rmpfr` package for multiple-precision evaluation of the series.

## Parameters

- `tau` is a fixed probability level in $(0,1)$.
- `mu_tau` is the corresponding quantile, positive.
- `sigma` is the positive parameter of the Akash base distribution.
- `alpha` is the positive power parameter obtained with
  `qeakash_alpha(mu_tau, sigma, tau)`,
  $\alpha = \log(\tau)/\log F(\mu_\tau;\sigma)$, where $F$ is the Akash
  distribution function.

The moment integral is

$$
\mu_r' = \alpha\,\frac{\sigma^3}{\sigma^2+2}\int_0^\infty x^r(1+x^2)e^{-\sigma x}F(x;\sigma)^{\alpha-1}\,dx .
$$

Writing the Akash survival function as $W(x)e^{-\sigma x}$, with
$W(x) = 1 + Ax + Bx^2$, $A = 2\sigma/(\sigma^2+2)$ and
$B = \sigma^2/(\sigma^2+2)$, the generalized binomial expansion of
$F(x;\sigma)^{\alpha-1}$ and term-by-term integration give the absolutely
convergent series

$$
\mu_r' = \frac{\alpha\sigma^3}{\sigma^2+2}\sum_{j=0}^{\infty}(-1)^j\binom{\alpha-1}{j}\sum_{m=0}^{2j} b_{j,m}\left[\frac{(r+m)!}{c_j^{\,r+m+1}} + \frac{(r+m+2)!}{c_j^{\,r+m+3}}\right],
\qquad c_j = (j+1)\sigma,
$$

where $b_{j,m}$ is the coefficient of $x^m$ in $W(x)^j$, obtained by the
recurrence $b_{0,0}=1$ and
$b_{j,m} = b_{j-1,m} + A\,b_{j-1,m-1} + B\,b_{j-1,m-2}$. When $\alpha = 1$,
that is, $\tau = F(\mu_\tau;\sigma)$, the distribution reduces to the Akash
distribution, with moments
$r!\,[\sigma^2+(r+1)(r+2)]/[\sigma^r(\sigma^2+2)]$.

## Example 1. Mean and variance

From the directory containing the script, load the functions into a
separate environment to avoid triggering the full comparison and file exports.

```r
qeakash <- new.env()
sys.source("QEAkash_Moments.R", envir = qeakash)

m1 <- qeakash$qeakash_moment_integral(r = 1, mu_tau = 2.5, sigma = 1.8, tau = 0.5)
m2 <- qeakash$qeakash_moment_integral(r = 2, mu_tau = 2.5, sigma = 1.8, tau = 0.5)
c(mean = m1, variance = m2 - m1^2)
```

The first raw moment is the mean, `2.6363`. The variance is the second raw
moment, `7.8917`, minus the square of the first.

## Example 2. Compare numerical methods

```r
cases <- data.frame(mu_tau = c(2.5, 0.5), sigma = c(1.8, 0.9))
comparison <- qeakash$qeakash_compare_moments(cases, taus = c(0.25, 0.75),
                                              moments = 1:2, J = 200)

comparison[c("tau", "mu_tau", "sigma", "alpha", "r", "moment_integral",
             "moment_series", "abs_difference")]
```

Each row represents a parameter pair, a probability level and a moment
order. `abs_difference` is $D_r = |\mu_r'^{\,I} - \mu_r'^{\,S}|$, the
absolute difference between quadrature and the truncated series.

| Function | Method |
|:--|:--|
| `qeakash_moment_integral(r, mu_tau, sigma, tau)` | Adaptive Gauss–Kronrod quadrature with `stats::integrate()`, split at the median. When $\alpha<1$, the substitution $v = x^\alpha$ removes the integrable singularity at the origin. This is the reference used in the comparisons. |
| `qeakash_moment_series(r, mu_tau, sigma, tau, J)` | Binomial series truncated at $j = J$, for one or more orders `r`. With `Rmpfr` installed, the sum is accumulated with 512 bits. Set `precision_bits = NULL` for the double-precision version. |
| `qeakash_moment_akash(r, sigma)` | Closed-form Akash moments, the case $\alpha = 1$. |
| `qeakash_compare_moments(cases, taus, moments, J)` | Both methods for each case, level and order, with the absolute differences. |

For $j<\alpha$, the binomial weights alternate in sign and grow with
$\alpha$. In double precision the series then loses digits through
cancellation: with $\mu_\tau=2.5$, $\sigma=1.8$ and $\tau=0.10$, the
relative error of the mean is about $10^{-8}$. Multiple precision removes
this source of error, so that $D_r$ measures truncation and quadrature
error. Truncation error is largest when $\alpha$ is small, where the series
converges more slowly.

## Run the full comparison

```sh
Rscript QEAkash_Moments.R
```

This evaluates moments of orders 1–4 for eight pairs
$(\mu_\tau, \sigma)$ = (0.5, 0.9), (1.5, 0.9), (2.5, 0.9), (0.5, 1.8),
(1.5, 1.8), (2.5, 1.8), (0.5, 2.9) and (1.5, 2.9), each crossed with
`tau = 0.10, 0.25, 0.50, 0.75, 0.90`, with the series truncated at
$J = 200$. The run takes several minutes. The script prints the 160 rows
and writes the following files to the current working directory:

- `qeakash_moment_comparison.csv` contains the estimates and discrepancies.
- `qeakash_moment_comparison.tex` contains a LaTeX table fragment. Load
  `booktabs` in the parent document.

Running the script again replaces these outputs. Scenario 6,
$\mu_\tau = 2.5$ and $\sigma = 1.8$, gives the following moments:

| `tau` | $\mu_1'$ | $\mu_2'$ | $\mu_3'$ | $\mu_4'$ | largest $D_r$ |
|:--|--:|--:|--:|--:|--:|
| 0.10 | 3.5474 | 13.4400 | 54.5104 | 237.2260 | 1.71e-13 |
| 0.25 | 3.1680 | 10.9265 | 41.0478 | 168.0442 | 8.53e-14 |
| 0.50 | 2.6363 | 7.8917 | 26.6181 | 100.5502 | 5.68e-14 |
| 0.75 | 1.9465 | 4.7783 | 14.1214 | 48.8426 | 6.20e-12 |
| 0.90 | 1.1950 | 2.3385 | 6.1205 | 19.7484 | 4.40e-07 |
