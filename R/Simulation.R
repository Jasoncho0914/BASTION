#' Simulate data with known calendar seasonality
#'
#' Generate observations, timestamps, and known components for calendar fitting
#' examples. See Details for the three simulation definitions.
#' @param scenario One of "daily", "biweekly", or "hourly"; default "daily".
#' @param seed Optional nonnegative integer random seed.
#' @details Each scenario generates y as the sum of its known components and
#'   independent zero-mean Gaussian noise. Timestamps use UTC.
#'
#'   "daily" generates daily observations from 2019 through 2021, including
#'   February 29, 2020. For t starting at zero, the trend is 50 + 0.005*t.
#'   With midpoint calendar phases w within the week and a within the year,
#'   weekly seasonality is 3*cos(2*pi*w), and yearly seasonality is
#'   8*sin(2*pi*a) + 2*cos(4*pi*a). Noise SD is 1.5.
#'
#'   "biweekly" generates 560 daily observations starting January 6, 2020.
#'   With t starting at zero, the trend is 100 + 0.01*t + 3*sin(2*pi*t/560),
#'   and seasonality is 6*sin(2*pi*t/14) + 2*cos(4*pi*t/14).
#'   Noise SD is 1.5.
#'
#'   "hourly" generates 2016 hourly observations starting January 1, 2024
#'   (12 full weeks). With t starting at zero, the trend is 100 + 0.002*t.
#'   With midpoint phase h within each day, daily seasonality is
#'   10*sin(2*pi*h) + 2*cos(4*pi*h). Monday-to-Sunday effects are
#'   c(2, 3, 3, 2, 1, -5, -6), constant within each day. Noise SD is 2.
#' @return A list with numeric observations `y`, POSIXct timestamps `dates`,
#'   a data frame `truth` containing the trend and seasonal components,
#'   the Gaussian `noise` vector, and its generating `noise_sd`.
#' @export
simulate_BASTION_calendar <- function(scenario=c("daily","biweekly","hourly"),seed=NULL) {
  scenario<-match.arg(scenario)
  if(!is.null(seed)) {
    if(!is.numeric(seed) || length(seed)!=1L || !is.finite(seed) ||
       seed<0 || seed>.Machine$integer.max || seed!=floor(seed))
      stop("seed must be a nonnegative integer or NULL")
    set.seed(as.integer(seed))
  }
  if(scenario=="daily") {
    dates<-seq(as.POSIXct("2019-01-01",tz="UTC"),
      as.POSIXct("2021-12-31",tz="UTC"),by="day")
    n<-length(dates);t<-0:(n-1L)
    year<-as.integer(format(dates,"%Y",tz="UTC"))
    year_start<-as.Date(paste0(year,"-01-01"))
    year_end<-as.Date(paste0(year+1L,"-01-01"))
    year_phase<-(as.numeric(as.Date(dates,tz="UTC")-year_start)+.5)/
      as.numeric(year_end-year_start)
    weekday<-as.integer(format(dates,"%u",tz="UTC"))-1L
    week_phase<-(weekday+.5)/7
    truth<-data.frame(trend=50+.005*t,
      seasonal_weekly=3*cos(2*pi*week_phase),
      seasonal_yearly=8*sin(2*pi*year_phase)+2*cos(4*pi*year_phase))
    noise_sd<-1.5
  } else if(scenario=="biweekly") {
    n<-14L*40L;t<-0:(n-1L)
    dates<-seq(as.POSIXct("2020-01-06",tz="UTC"),by="day",length.out=n)
    truth<-data.frame(trend=100+.01*t+3*sin(2*pi*t/n),
      seasonal_biweekly=6*sin(2*pi*t/14)+2*cos(4*pi*t/14))
    noise_sd<-1.5
  } else {
    n<-24L*7L*12L;t<-0:(n-1L)
    dates<-seq(as.POSIXct("2024-01-01",tz="UTC"),by="hour",length.out=n)
    hour_phase<-((t%%24)+.5)/24
    weekday<-(t%/%24)%%7
    weekday_effect<-c(2,3,3,2,1,-5,-6)
    truth<-data.frame(trend=100+.002*t,
      seasonal_daily=10*sin(2*pi*hour_phase)+2*cos(4*pi*hour_phase),
      seasonal_weekly=weekday_effect[weekday+1L])
    noise_sd<-2
  }
  noise<-stats::rnorm(n,sd=noise_sd)
  list(y=rowSums(truth)+noise,dates=dates,truth=truth,noise=noise,noise_sd=noise_sd)
}

#' Simulate data with known fixed seasonal periods
#'
#' Generate observations and known components for fixed-period fitting examples.
#' @param scenario One of "daily", "biweekly", or "hourly"; default "daily".
#' @param seed Optional nonnegative integer random seed.
#' @details Each scenario generates y as the sum of its known components and
#'   independent zero-mean Gaussian noise. Time is indexed by observation number.
#'
#'   "daily" generates 1095 observations with fixed periods 7 and 365.
#'   With t starting at zero, the trend is 50 + 0.005*t. Weekly seasonality is
#'   3*cos(2*pi*(t+0.5)/7), and yearly seasonality is
#'   8*sin(2*pi*(t+0.5)/365) + 2*cos(4*pi*(t+0.5)/365). Noise SD is 1.5.
#'
#'   "biweekly" generates 560 observations with period 14. With t starting at
#'   zero, the trend is 100 + 0.01*t + 3*sin(2*pi*t/560), and seasonality is
#'   6*sin(2*pi*t/14) + 2*cos(4*pi*t/14). Noise SD is 1.5.
#'
#'   "hourly" generates 2016 observations with periods 24 and 168. With t
#'   starting at zero, the trend is 100 + 0.002*t. Using midpoint phase h within
#'   each 24-observation cycle, daily seasonality is
#'   10*sin(2*pi*h) + 2*cos(4*pi*h). Weekly effects are c(2, 3, 3, 2, 1, -5, -6),
#'   each held for 24 observations and repeated every 168. Noise SD is 2.
#' @return A list with numeric observations `y`, observation indices `time`,
#'   a data frame `truth` containing the trend and seasonal components,
#'   the Gaussian `noise` vector, and its generating `noise_sd`.
#' @seealso [fit_BASTION()], [simulate_BASTION_calendar()]
#' @export
simulate_BASTION <- function(scenario=c("daily","biweekly","hourly"),seed=NULL) {
  scenario<-match.arg(scenario)
  if(!is.null(seed)) {
    if(!is.numeric(seed) || length(seed)!=1L || !is.finite(seed) ||
       seed<0 || seed>.Machine$integer.max || seed!=floor(seed))
      stop("seed must be a nonnegative integer or NULL")
    set.seed(as.integer(seed))
  }
  if(scenario=="daily") {
    n<-3L*365L;t<-0:(n-1L)
    truth<-data.frame(trend=50+.005*t,
      seasonal_7=3*cos(2*pi*(t+.5)/7),
      seasonal_365=8*sin(2*pi*(t+.5)/365)+2*cos(4*pi*(t+.5)/365))
    noise_sd<-1.5
  } else if(scenario=="biweekly") {
    n<-14L*40L;t<-0:(n-1L)
    truth<-data.frame(trend=100+.01*t+3*sin(2*pi*t/n),
      seasonal_14=6*sin(2*pi*t/14)+2*cos(4*pi*t/14))
    noise_sd<-1.5
  } else {
    n<-24L*7L*12L;t<-0:(n-1L)
    hour_phase<-((t%%24)+.5)/24
    weekday<-(t%/%24)%%7
    weekday_effect<-c(2,3,3,2,1,-5,-6)
    truth<-data.frame(trend=100+.002*t,
      seasonal_24=10*sin(2*pi*hour_phase)+2*cos(4*pi*hour_phase),
      seasonal_168=weekday_effect[weekday+1L])
    noise_sd<-2
  }
  noise<-stats::rnorm(n,sd=noise_sd)
  list(y=rowSums(truth)+noise,time=seq_len(n),truth=truth,noise=noise,noise_sd=noise_sd)
}
