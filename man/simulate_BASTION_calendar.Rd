\name{simulate_BASTION_calendar}
\alias{simulate_BASTION_calendar}
\title{Simulate data with known calendar seasonality}
\description{Generate observations, timestamps, and known components for calendar fitting examples.}
\usage{simulate_BASTION_calendar(scenario = c("daily", "biweekly", "hourly"), seed = NULL)}
\arguments{
\item{scenario}{One of "daily", "biweekly", or "hourly"; default "daily".}
\item{seed}{Optional nonnegative integer random seed.}
}
\details{
Each scenario generates y as the sum of its known components and independent zero-mean Gaussian noise. Timestamps use UTC.

\code{"daily"} generates daily observations from 2019 through 2021, including February 29, 2020. With t starting at zero, the trend is 50 + 0.005*t. Using midpoint calendar phases w within the week and a within the year, weekly seasonality is 3*cos(2*pi*w), and yearly seasonality is 8*sin(2*pi*a) + 2*cos(4*pi*a). Noise SD is 1.5.

\code{"biweekly"} generates 560 daily observations starting January 6, 2020. With t starting at zero, the trend is 100 + 0.01*t + 3*sin(2*pi*t/560), and seasonality is 6*sin(2*pi*t/14) + 2*cos(4*pi*t/14). Noise SD is 1.5.

\code{"hourly"} generates 2016 hourly observations starting January 1, 2024 (12 full weeks). With t starting at zero, the trend is 100 + 0.002*t. Using midpoint phase h within each day, daily seasonality is 10*sin(2*pi*h) + 2*cos(4*pi*h). Monday-to-Sunday effects are c(2, 3, 3, 2, 1, -5, -6), constant within each day. Noise SD is 2.
}
\value{A list with numeric observations \code{y}, POSIXct timestamps \code{dates}, a data frame \code{truth} containing the trend and seasonal components, the Gaussian \code{noise} vector, and its generating \code{noise_sd}.}
\seealso{\code{\link{fit_BASTION_calendar}}}
