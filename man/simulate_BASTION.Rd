\name{simulate_BASTION}
\alias{simulate_BASTION}
\title{Simulate data with known fixed seasonal periods}
\description{Generate observations and known components for fixed-period fitting examples.}
\usage{simulate_BASTION(scenario = c("daily", "biweekly", "hourly"), seed = NULL)}
\arguments{
\item{scenario}{One of "daily", "biweekly", or "hourly"; default "daily".}
\item{seed}{Optional nonnegative integer random seed.}
}
\details{
Each scenario generates y as the sum of its known components and independent zero-mean Gaussian noise. Time is indexed by observation number.

\code{"daily"} generates 1095 observations with fixed periods 7 and 365. With t starting at zero, the trend is 50 + 0.005*t. Weekly seasonality is 3*cos(2*pi*(t+0.5)/7), and yearly seasonality is 8*sin(2*pi*(t+0.5)/365) + 2*cos(4*pi*(t+0.5)/365). Noise SD is 1.5.

\code{"biweekly"} generates 560 observations with period 14. With t starting at zero, the trend is 100 + 0.01*t + 3*sin(2*pi*t/560), and seasonality is 6*sin(2*pi*t/14) + 2*cos(4*pi*t/14). Noise SD is 1.5.

\code{"hourly"} generates 2016 observations with periods 24 and 168. With t starting at zero, the trend is 100 + 0.002*t. Using midpoint phase h within each 24-observation cycle, daily seasonality is 10*sin(2*pi*h) + 2*cos(4*pi*h). Weekly effects are c(2, 3, 3, 2, 1, -5, -6), each held for 24 observations and repeated every 168. Noise SD is 2.
}
\value{A list with numeric observations \code{y}, observation indices \code{time}, a data frame \code{truth} containing the trend and seasonal components, the Gaussian \code{noise} vector, and its generating \code{noise_sd}.}
\seealso{\code{\link{fit_BASTION}}, \code{\link{simulate_BASTION_calendar}}}
