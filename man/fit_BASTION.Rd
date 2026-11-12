\name{fit_BASTION}
\alias{fit_BASTION}
\title{BASTION decomposition}
\description{Bayesian trend and seasonal decomposition with locally adaptive priors.}
\usage{
fit_BASTION(y, Ks, X = NULL, Outlier = FALSE, cl = 0.95,
  sparse = FALSE, obsSV = "const", nchains = 2, nsave = 1000,
  nburn = 1000, nskip = 4, verbose = TRUE, save_samples = FALSE,
  D_trend = 2L)
}
\arguments{
\item{y}{Numeric observation vector or univariate time series.}
\item{Ks}{List of fixed seasonal periods.}
\item{X}{Optional regression design matrix.}
\item{Outlier}{Include outliers with a horseshoe+ prior at all observations (default FALSE).}
\item{cl}{Pointwise credible interval probability.}
\item{sparse}{Additional trend shrinkage option.}
\item{obsSV}{Observation variance model: "const" or "SV".}
\item{nchains}{Number of chains.}
\item{nsave}{Saved draws per chain.}
\item{nburn}{Burn-in sweeps per chain.}
\item{nskip}{Discarded sweeps between saves; save every nskip+1 sweeps.}
\item{verbose}{Display progress.}
\item{save_samples}{Include posterior draws under \code{draws}.}
\item{D_trend}{Trend difference order, 1 or 2 (default 2).}
}
\value{
An object of class \code{bastion_fit}. Let T be the number of observations,
M the number of decomposition columns (including residual), D the total saved draws,
K a seasonal period, and p the number of regression coefficients.
\describe{
\item{time}{Observation times, length T. Uses \code{time(y)} for a time series and observation numbers otherwise.}
\item{observed}{Supplied observations, length T.}
\item{components}{T by M data frame of posterior means. Columns are \code{trend}, \code{seasonal_<period>}, \code{regression}/\code{outlier} when fitted, and \code{residual}.}
\item{intervals}{\code{lower} and \code{upper}, T by M data frames with matching component names and order.}
\item{seasonality}{One named entry per seasonal component. Each contains \code{recurring}, a K by 3 data frame with mean, lower, upper; and \code{deviation}, a T by 3 data frame with mean, lower, upper.}
\item{regression}{When fitted, \code{coefficients} is a p by 4 data frame with term, mean, lower, upper. Unnamed covariates use beta_1, beta_2, etc.}
\item{noise}{\code{sd} has T rows for SV or one row for constant SD, with mean, lower, upper. Under SV, \code{parameters} has three rows with parameter, mean, lower, upper.}
\item{metadata}{Fitting metadata includes \code{cl}, the credible-interval probability, and \code{alpha = 1 - cl}. Input seasonal periods, fitting arguments, and recurring positions under mapping$recurring_positions.}
\item{draws}{Present only when \code{save_samples=TRUE}. \code{index} is a D by 3 data frame (draw, chain, iteration); \code{components} is D by T by M, including residuals. Each \code{seasonality} entry contains D by K \code{recurring} and D by T \code{deviation} matrices. Optional \code{regression$coefficients} is D by p. \code{noise$sd} is D by T under SV or D by 1 otherwise; SV \code{noise$parameters} is D by 3 (mu, phi, innovation_sd).}
}
}
\details{
BASTION decomposes a time series into a locally adaptive trend, one or more seasonal components, and optional regression and outlier components. The components are estimated jointly using Bayesian inference.

Supply seasonal periods in numbers of observations through \code{Ks}. For example, use \code{list(12)} for yearly seasonality in monthly data, or \code{list(7, 365)} for weekly and yearly seasonality in daily data. Each seasonal component combines a common recurring pattern with deviations that allow the pattern to evolve across cycles.

The trend uses second differences by default (\code{D_trend = 2}). Use \code{D_trend = 1} for first differences. Setting \code{sparse = TRUE} applies stronger shrinkage to trend differences while retaining local adaptivity.

Use \code{obsSV = "const"} for constant observation variance or \code{obsSV = "SV"} to estimate observation variability over time. Set \code{Outlier = TRUE} to estimate a separate outlier component, and supply \code{X} to include regression effects.

Component estimates and observation standard deviations are returned in the units of the supplied data. Posterior means and credible intervals at the level specified by \code{cl} are returned whether or not posterior draws are retained. Set \code{save_samples = TRUE} to retain the draws.

Use \code{print()} for an overview of the fit, \code{summary()} for the component estimates, \code{plot()} for the decomposition, and \code{plot_component()} to inspect an individual component with credible intervals.
}

\seealso{\code{\link{simulate_BASTION}}}
\examples{
\dontrun{
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
}
}
