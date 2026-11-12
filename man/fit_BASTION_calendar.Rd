\name{fit_BASTION_calendar}
\alias{fit_BASTION_calendar}
\title{BASTION decomposition with calendar seasonality}
\description{Bayesian trend and seasonal decomposition with locally adaptive priors.}
\usage{
fit_BASTION_calendar(y, dates, data_unit, seasonalities, X = NULL,
  Outlier = FALSE, cl = .95, sparse = FALSE, obsSV = "const",
  nchains = 2L, nsave = 1000L, nburn = 1000L, nskip = 4L,
  verbose = TRUE, save_samples = FALSE, D_trend = 2L,
  tz = NULL, tolerance_seconds = NULL, deviation_periods = NULL)
}
\arguments{
\item{y}{Numeric observations; NA allowed.}
\item{dates}{POSIXct timestamps matching the observations.}
\item{data_unit}{Collection interval list with unit, optional every and anchor.}
\item{seasonalities}{Named seasonal specification lists with unit, optional every, anchor and distance_scale.}
\item{X}{Optional regression design matrix.}
\item{Outlier}{Include outliers with a horseshoe+ prior at all observations (default FALSE).}
\item{cl}{Pointwise equal-tailed credible interval probability.}
\item{sparse}{Additional trend shrinkage option; not a matrix storage switch.}
\item{obsSV}{Observation variance model: "const" or "SV".}
\item{nchains}{Number of chains.}
\item{nsave}{Saved draws per chain.}
\item{nburn}{Burn-in sweeps per chain.}
\item{nskip}{Discarded sweeps between saves; save every nskip+1 sweeps.}
\item{verbose}{Display progress.}
\item{save_samples}{Return posterior draws.}
\item{D_trend}{Trend difference order, 1 or 2 (default 2).}
\item{tz}{Optional explicit time zone.}
\item{tolerance_seconds}{Timestamp matching tolerance passed to preprocessing.}
\item{deviation_periods}{Optional named fixed-lag overrides. Supply these for modal ties or no complete cycles.}
}
\details{
BASTION decomposes a time series into a locally adaptive trend, one or more seasonal components, and optional regression and outlier components. The calendar fitter uses observation dates to align seasonal patterns with calendar cycles. All components are estimated jointly using Bayesian inference.

Supply observations through \code{y} and matching \code{POSIXct} timestamps through \code{dates}. Use \code{data_unit} to describe the observation interval, such as \code{list(unit = "day")} for daily data or \code{list(unit = "hour")} for hourly data.

Specify seasonal patterns as a named list through \code{seasonalities}. For example, \code{list(weekly = list(unit = "week"), yearly = list(unit = "year"))} fits weekly and yearly seasonality together. Calendar alignment accommodates different cycle lengths, such as months with different numbers of days and leap years. Each seasonal component combines a smooth recurring pattern with deviations that allow the pattern to evolve across cycles.

The trend uses second differences by default (\code{D_trend = 2}). Use \code{D_trend = 1} for first differences. Setting \code{sparse = TRUE} applies stronger shrinkage to trend differences while retaining local adaptivity.

Use \code{obsSV = "const"} for constant observation variance or \code{obsSV = "SV"} to estimate observation variability over time. Set \code{Outlier = TRUE} to estimate a separate outlier component, and supply \code{X} to include regression effects.

Missing observations are estimated during fitting. The returned timeline includes all scheduled timestamps, and \code{imputation} identifies those with missing observations. Component estimates and observation standard deviations are returned in the units of the supplied data. Posterior means and credible intervals at the level specified by \code{cl} are returned whether or not posterior draws are retained. Set \code{save_samples = TRUE} to retain the draws.

Use \code{print()} for an overview of the fit, \code{summary()} for the component estimates, \code{plot()} for the decomposition, and \code{plot_component()} to inspect an individual component with credible intervals. Seasonal component names follow the names supplied in \code{seasonalities}, such as \code{seasonal_weekly} and \code{seasonal_yearly}.
}
\value{
An object of class \code{bastion_fit}, with the same component, interval, seasonality, regression, noise and optional draw structure as \code{fit_BASTION}. Let T be all scheduled timestamps, D saved draws, M decomposition columns (including residual), J native phase-grid positions, and H missing collection slots.
\describe{
\item{time, observed}{All scheduled timestamps in chronological order with the resolved time zone, and observations, length T. Missing responses are NA, including omitted timestamps.}
\item{components}{T by M posterior-mean data frame, with trend, seasonal_<supplied name>, regression/outlier when fitted, and residual.}
\item{intervals}{T by M lower and upper data frames with matching component names/order.}
\item{seasonality}{Each named entry has recurring, a J by 3 data frame (mean, lower, upper), and deviation, a T by 3 data frame (mean, lower, upper). Native phase coordinates in [0,1) are stored in metadata$mapping$recurring_positions.}
\item{regression}{Optional coefficient summaries, as in fit_BASTION.}
\item{noise}{SD summaries and SV parameter summaries when applicable, as in fit_BASTION. Residuals are included in components and intervals.}
\item{metadata}{Fitting metadata includes \code{cl}, the credible-interval probability, and \code{alpha = 1 - cl}. Input specifications (including collection unit/anchor and time zone), fitting arguments, and compact seasonality mappings under mapping$seasonality. Each seasonality has its unit, every, anchor, phase_index (length T), pattern_id (length T), and patterns (a list of full-cycle grid-index vectors). These indices refer to rows of its recurring table.}
\item{imputation}{POSIXct vector of H timestamps with missing responses, including explicit missing responses and omitted collection timestamps. All component estimates and intervals remain on the full timeline.}
\item{draws}{Only when save_samples=TRUE. Components are D by T by M; recurring patterns D by J; deviations are D by T, and residuals are included in the component array. Noise SD is D by T under SV or D by 1 otherwise; SV parameters are D by 3. Regression coefficients and draw/chain/iteration indices follow fit_BASTION. Completed-response draws can be recovered by summing components, including residual, within each draw.}
}
}

\section{Supplying observation dates and intervals}{
\code{y} must be a numeric vector, and \code{dates} must be a \code{POSIXct} vector of the same length, with a valid time zone. Each timestamp identifies the observation in the corresponding position of \code{y}. Convert character dates before calling the fitter, for example with \code{as.POSIXct(date_text, tz = "UTC")}.

\code{data_unit} describes how often observations are collected. Supply a list containing \code{unit} and, optionally, \code{every} and \code{anchor}. For example, \code{list(unit = "day")} means one observation per day, \code{list(unit = "day", every = 2)} means one observation every two days, and \code{list(unit = "minute", every = 15)} means one observation every 15 minutes. The default value of \code{every} is 1.

Use the intended observation interval even when some observations are missing. For daily data with missing days, keep \code{data_unit = list(unit = "day")}. Missing responses may be supplied as \code{NA}; omitted timestamps are also recognized as missing positions in that interval. If \code{X} is supplied, include every timestamp and its covariates, using \code{NA} in \code{y} for missing responses.

Allowed units are \code{"millisecond"}, \code{"second"}, \code{"minute"}, \code{"hour"}, \code{"day"}, \code{"week"}, \code{"month"}, \code{"quarter"}, and \code{"year"}. Day, week, month, quarter, and year use calendar intervals and require a positive integer \code{every}. Subdaily units use elapsed time and allow positive fractional values, such as \code{list(unit = "second", every = 0.01)} for observations every 10 milliseconds.
}

\section{Specifying seasonal patterns}{
\code{seasonalities} is a named list with one list per seasonal component. Each component has \code{unit}, optional \code{every} (default 1), optional \code{anchor}, and optional \code{distance_scale} (default \code{TRUE}). The allowed units are the same as for \code{data_unit}. The list names label the results; for example, \code{biweekly = list(unit = "week", every = 2)} creates a component named \code{seasonal_biweekly}.

Here \code{unit} and \code{every} describe how often a seasonal pattern repeats. With daily observations, \code{list(unit = "week")} repeats every week, \code{list(unit = "week", every = 2)} repeats every two weeks, and \code{list(unit = "month")} repeats every calendar month. With hourly observations, \code{list(unit = "day")} repeats every day. Supplying multiple entries estimates those seasonal components together.

Keep the observation interval and the seasonal period separate: daily observations with a biweekly pattern use \code{data_unit = list(unit = "day")} and \code{seasonalities = list(biweekly = list(unit = "week", every = 2))}. The seasonal cycle must contain enough observations to estimate its shape; a daily pattern cannot be estimated from one observation per day.

\code{anchor} sets a reference boundary for an interval or seasonal cycle and must be a single \code{POSIXct} timestamp. The observation interval defaults to the earliest supplied timestamp. Seasonal boundaries default to midnight on January 1, 1970, except weekly cycles use Monday, January 5, 1970, in the selected time zone. For a two-week cycle, set a custom anchor to specify which week starts the cycle.

\code{distance_scale = TRUE} accounts for spacing between seasonal positions when smoothing a grid-based recurring pattern. Set it to \code{FALSE} to use unscaled second differences for that pattern. This option is supplied separately within each seasonal specification; it does not change the observation interval or the number of days in a calendar cycle.
}
\seealso{\code{\link{simulate_BASTION_calendar}}}
\examples{
\dontrun{
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
}
}
