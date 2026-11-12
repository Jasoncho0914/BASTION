# BASTION decomposition

Bayesian trend and seasonal decomposition with locally adaptive priors.

## Usage

``` r
fit_BASTION(y, Ks, X = NULL, Outlier = FALSE, cl = 0.95,
  sparse = FALSE, obsSV = "const", nchains = 2, nsave = 1000,
  nburn = 1000, nskip = 4, verbose = TRUE, save_samples = FALSE,
  D_trend = 2L)
```

## Arguments

- y:

  Numeric observation vector or univariate time series.

- Ks:

  List of fixed seasonal periods.

- X:

  Optional regression design matrix.

- Outlier:

  Include outliers with a horseshoe+ prior at all observations (default
  FALSE).

- cl:

  Pointwise credible interval probability.

- sparse:

  Additional trend shrinkage option.

- obsSV:

  Observation variance model: "const" or "SV".

- nchains:

  Number of chains.

- nsave:

  Saved draws per chain.

- nburn:

  Burn-in sweeps per chain.

- nskip:

  Discarded sweeps between saves; save every nskip+1 sweeps.

- verbose:

  Display progress.

- save_samples:

  Include posterior draws under `draws`.

- D_trend:

  Trend difference order, 1 or 2 (default 2).

## Value

An object of class `bastion_fit`. Let T be the number of observations, M
the number of decomposition columns (including residual), D the total
saved draws, K a seasonal period, and p the number of regression
coefficients.

- time:

  Observation times, length T. Uses `time(y)` for a time series and
  observation numbers otherwise.

- observed:

  Supplied observations, length T.

- components:

  T by M data frame of posterior means. Columns are `trend`,
  `seasonal_<period>`, `regression`/`outlier` when fitted, and
  `residual`.

- intervals:

  `lower` and `upper`, T by M data frames with matching component names
  and order.

- seasonality:

  One named entry per seasonal component. Each contains `recurring`, a K
  by 3 data frame with mean, lower, upper; and `deviation`, a T by 3
  data frame with mean, lower, upper.

- regression:

  When fitted, `coefficients` is a p by 4 data frame with term, mean,
  lower, upper. Unnamed covariates use beta_1, beta_2, etc.

- noise:

  `sd` has T rows for SV or one row for constant SD, with mean, lower,
  upper. Under SV, `parameters` has three rows with parameter, mean,
  lower, upper.

- metadata:

  Fitting metadata includes `cl`, the credible-interval probability, and
  `alpha = 1 - cl`. Input seasonal periods, fitting arguments, and
  recurring positions under mapping\$recurring_positions.

- draws:

  Present only when `save_samples=TRUE`. `index` is a D by 3 data frame
  (draw, chain, iteration); `components` is D by T by M, including
  residuals. Each `seasonality` entry contains D by K `recurring` and D
  by T `deviation` matrices. Optional `regression$coefficients` is D
  by p. `noise$sd` is D by T under SV or D by 1 otherwise; SV
  `noise$parameters` is D by 3 (mu, phi, innovation_sd).

## Details

BASTION decomposes a time series into a locally adaptive trend, one or
more seasonal components, and optional regression and outlier
components. The components are estimated jointly using Bayesian
inference.

Supply seasonal periods in numbers of observations through `Ks`. For
example, use `list(12)` for yearly seasonality in monthly data, or
`list(7, 365)` for weekly and yearly seasonality in daily data. Each
seasonal component combines a common recurring pattern with deviations
that allow the pattern to evolve across cycles.

The trend uses second differences by default (`D_trend = 2`). Use
`D_trend = 1` for first differences. Setting `sparse = TRUE` applies
stronger shrinkage to trend differences while retaining local
adaptivity.

Use `obsSV = "const"` for constant observation variance or
`obsSV = "SV"` to estimate observation variability over time. Set
`Outlier = TRUE` to estimate a separate outlier component, and supply
`X` to include regression effects.

Component estimates and observation standard deviations are returned in
the units of the supplied data. Posterior means and credible intervals
at the level specified by `cl` are returned whether or not posterior
draws are retained. Set `save_samples = TRUE` to retain the draws.

Use [`print()`](https://rdrr.io/r/base/print.html) for an overview of
the fit, [`summary()`](https://rdrr.io/r/base/summary.html) for the
component estimates,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for the
decomposition, and
[`plot_component()`](https://jasoncho0914.github.io/BASTION/reference/plot_component.md)
to inspect an individual component with credible intervals.

## See also

[`simulate_BASTION`](https://jasoncho0914.github.io/BASTION/reference/simulate_BASTION.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Daily observations with fixed weekly and yearly periods
sim <- simulate_BASTION("daily", seed = 20)
fit_daily <- fit_BASTION(
  y = sim$y,
  Ks = list(7, 365),
  obsSV = "const", nchains = 1
)
plot(fit_daily)
plot_component(fit_daily, "seasonal_365", part = "recurring")

# A pattern repeating every 14 observations
sim <- simulate_BASTION("biweekly", seed = 21)
fit_biweekly <- fit_BASTION(
  y = sim$y,
  Ks = list(14),
  obsSV = "const", nchains = 1
)
plot_component(fit_biweekly, "seasonal_14")

# Hourly observations with daily and weekly periods
sim <- simulate_BASTION("hourly", seed = 22)
fit_hourly <- fit_BASTION(
  y = sim$y,
  Ks = list(24, 168),
  obsSV = "const", nchains = 1
)
plot(fit_hourly)
} # }
```
