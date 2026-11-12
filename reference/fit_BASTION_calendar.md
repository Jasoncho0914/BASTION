# BASTION decomposition with calendar seasonality

Bayesian trend and seasonal decomposition with locally adaptive priors.

## Usage

``` r
fit_BASTION_calendar(y, dates, data_unit, seasonalities, X = NULL,
  Outlier = FALSE, cl = .95, sparse = FALSE, obsSV = "const",
  nchains = 2L, nsave = 1000L, nburn = 1000L, nskip = 4L,
  verbose = TRUE, save_samples = FALSE, D_trend = 2L,
  tz = NULL, tolerance_seconds = NULL, deviation_periods = NULL)
```

## Arguments

- y:

  Numeric observations; NA allowed.

- dates:

  POSIXct timestamps matching the observations.

- data_unit:

  Collection interval list with unit, optional every and anchor.

- seasonalities:

  Named seasonal specification lists with unit, optional every, anchor
  and distance_scale.

- X:

  Optional regression design matrix.

- Outlier:

  Include outliers with a horseshoe+ prior at all observations (default
  FALSE).

- cl:

  Pointwise equal-tailed credible interval probability.

- sparse:

  Additional trend shrinkage option; not a matrix storage switch.

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

  Return posterior draws.

- D_trend:

  Trend difference order, 1 or 2 (default 2).

- tz:

  Optional explicit time zone.

- tolerance_seconds:

  Timestamp matching tolerance passed to preprocessing.

- deviation_periods:

  Optional named fixed-lag overrides. Supply these for modal ties or no
  complete cycles.

## Details

BASTION decomposes a time series into a locally adaptive trend, one or
more seasonal components, and optional regression and outlier
components. The calendar fitter uses observation dates to align seasonal
patterns with calendar cycles. All components are estimated jointly
using Bayesian inference.

Supply observations through `y` and matching `POSIXct` timestamps
through `dates`. Use `data_unit` to describe the observation interval,
such as `list(unit = "day")` for daily data or `list(unit = "hour")` for
hourly data.

Specify seasonal patterns as a named list through `seasonalities`. For
example,
`list(weekly = list(unit = "week"), yearly = list(unit = "year"))` fits
weekly and yearly seasonality together. Calendar alignment accommodates
different cycle lengths, such as months with different numbers of days
and leap years. Each seasonal component combines a smooth recurring
pattern with deviations that allow the pattern to evolve across cycles.

The trend uses second differences by default (`D_trend = 2`). Use
`D_trend = 1` for first differences. Setting `sparse = TRUE` applies
stronger shrinkage to trend differences while retaining local
adaptivity.

Use `obsSV = "const"` for constant observation variance or
`obsSV = "SV"` to estimate observation variability over time. Set
`Outlier = TRUE` to estimate a separate outlier component, and supply
`X` to include regression effects.

Missing observations are estimated during fitting. The returned timeline
includes all scheduled timestamps, and `imputation` identifies those
with missing observations. Component estimates and observation standard
deviations are returned in the units of the supplied data. Posterior
means and credible intervals at the level specified by `cl` are returned
whether or not posterior draws are retained. Set `save_samples = TRUE`
to retain the draws.

Use [`print()`](https://rdrr.io/r/base/print.html) for an overview of
the fit, [`summary()`](https://rdrr.io/r/base/summary.html) for the
component estimates,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for the
decomposition, and
[`plot_component()`](https://jasoncho0914.github.io/BASTION/reference/plot_component.md)
to inspect an individual component with credible intervals. Seasonal
component names follow the names supplied in `seasonalities`, such as
`seasonal_weekly` and `seasonal_yearly`.

## Value

An object of class `bastion_fit`, with the same component, interval,
seasonality, regression, noise and optional draw structure as
`fit_BASTION`. Let T be all scheduled timestamps, D saved draws, M
decomposition columns (including residual), J native phase-grid
positions, and H missing collection slots.

- time, observed:

  All scheduled timestamps in chronological order with the resolved time
  zone, and observations, length T. Missing responses are NA, including
  omitted timestamps.

- components:

  T by M posterior-mean data frame, with trend, seasonal\_\<supplied
  name\>, regression/outlier when fitted, and residual.

- intervals:

  T by M lower and upper data frames with matching component
  names/order.

- seasonality:

  Each named entry has recurring, a J by 3 data frame (mean, lower,
  upper), and deviation, a T by 3 data frame (mean, lower, upper).
  Native phase coordinates in \[0,1) are stored in
  metadata\$mapping\$recurring_positions.

- regression:

  Optional coefficient summaries, as in fit_BASTION.

- noise:

  SD summaries and SV parameter summaries when applicable, as in
  fit_BASTION. Residuals are included in components and intervals.

- metadata:

  Fitting metadata includes `cl`, the credible-interval probability, and
  `alpha = 1 - cl`. Input specifications (including collection
  unit/anchor and time zone), fitting arguments, and compact seasonality
  mappings under mapping\$seasonality. Each seasonality has its unit,
  every, anchor, phase_index (length T), pattern_id (length T), and
  patterns (a list of full-cycle grid-index vectors). These indices
  refer to rows of its recurring table.

- imputation:

  POSIXct vector of H timestamps with missing responses, including
  explicit missing responses and omitted collection timestamps. All
  component estimates and intervals remain on the full timeline.

- draws:

  Only when save_samples=TRUE. Components are D by T by M; recurring
  patterns D by J; deviations are D by T, and residuals are included in
  the component array. Noise SD is D by T under SV or D by 1 otherwise;
  SV parameters are D by 3. Regression coefficients and
  draw/chain/iteration indices follow fit_BASTION. Completed-response
  draws can be recovered by summing components, including residual,
  within each draw.

## Supplying observation dates and intervals

`y` must be a numeric vector, and `dates` must be a `POSIXct` vector of
the same length, with a valid time zone. Each timestamp identifies the
observation in the corresponding position of `y`. Convert character
dates before calling the fitter, for example with
`as.POSIXct(date_text, tz = "UTC")`.

`data_unit` describes how often observations are collected. Supply a
list containing `unit` and, optionally, `every` and `anchor`. For
example, `list(unit = "day")` means one observation per day,
`list(unit = "day", every = 2)` means one observation every two days,
and `list(unit = "minute", every = 15)` means one observation every 15
minutes. The default value of `every` is 1.

Use the intended observation interval even when some observations are
missing. For daily data with missing days, keep
`data_unit = list(unit = "day")`. Missing responses may be supplied as
`NA`; omitted timestamps are also recognized as missing positions in
that interval. If `X` is supplied, include every timestamp and its
covariates, using `NA` in `y` for missing responses.

Allowed units are `"millisecond"`, `"second"`, `"minute"`, `"hour"`,
`"day"`, `"week"`, `"month"`, `"quarter"`, and `"year"`. Day, week,
month, quarter, and year use calendar intervals and require a positive
integer `every`. Subdaily units use elapsed time and allow positive
fractional values, such as `list(unit = "second", every = 0.01)` for
observations every 10 milliseconds.

## Specifying seasonal patterns

`seasonalities` is a named list with one list per seasonal component.
Each component has `unit`, optional `every` (default 1), optional
`anchor`, and optional `distance_scale` (default `TRUE`). The allowed
units are the same as for `data_unit`. The list names label the results;
for example, `biweekly = list(unit = "week", every = 2)` creates a
component named `seasonal_biweekly`.

Here `unit` and `every` describe how often a seasonal pattern repeats.
With daily observations, `list(unit = "week")` repeats every week,
`list(unit = "week", every = 2)` repeats every two weeks, and
`list(unit = "month")` repeats every calendar month. With hourly
observations, `list(unit = "day")` repeats every day. Supplying multiple
entries estimates those seasonal components together.

Keep the observation interval and the seasonal period separate: daily
observations with a biweekly pattern use
`data_unit = list(unit = "day")` and
`seasonalities = list(biweekly = list(unit = "week", every = 2))`. The
seasonal cycle must contain enough observations to estimate its shape; a
daily pattern cannot be estimated from one observation per day.

`anchor` sets a reference boundary for an interval or seasonal cycle and
must be a single `POSIXct` timestamp. The observation interval defaults
to the earliest supplied timestamp. Seasonal boundaries default to
midnight on January 1, 1970, except weekly cycles use Monday, January 5,
1970, in the selected time zone. For a two-week cycle, set a custom
anchor to specify which week starts the cycle.

`distance_scale = TRUE` accounts for spacing between seasonal positions
when smoothing a grid-based recurring pattern. Set it to `FALSE` to use
unscaled second differences for that pattern. This option is supplied
separately within each seasonal specification; it does not change the
observation interval or the number of days in a calendar cycle.

## See also

[`simulate_BASTION_calendar`](https://jasoncho0914.github.io/BASTION/reference/simulate_BASTION_calendar.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Daily observations with weekly and yearly patterns
sim <- simulate_BASTION_calendar("daily", seed = 20)
fit_daily <- fit_BASTION_calendar(
  y = sim$y,
  dates = sim$dates,
  data_unit = list(unit = "day"),
  seasonalities = list(
    weekly = list(unit = "week"),
    yearly = list(unit = "year")
  ),
  obsSV = "const", nchains = 1
)
plot(fit_daily)

# Daily observations; the pattern repeats every two weeks
sim <- simulate_BASTION_calendar("biweekly", seed = 21)
fit_biweekly <- fit_BASTION_calendar(
  y = sim$y,
  dates = sim$dates,
  data_unit = list(unit = "day"),
  seasonalities = list(
    biweekly = list(unit = "week", every = 2, anchor = sim$dates[1])
  ),
  obsSV = "const", nchains = 1
)
plot_component(fit_biweekly, "seasonal_biweekly")

# Hourly observations with daily and weekly patterns
sim <- simulate_BASTION_calendar("hourly", seed = 22)
fit_hourly <- fit_BASTION_calendar(
  y = sim$y,
  dates = sim$dates,
  data_unit = list(unit = "hour"),
  seasonalities = list(
    daily = list(unit = "day"),
    weekly = list(unit = "week")
  ),
  obsSV = "const", nchains = 1
)
plot(fit_hourly)
} # }
```
