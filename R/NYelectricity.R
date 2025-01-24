#' Daily Average Electricity Demand
#'
#' A dataset containing daily average electricity demand in New York State, aggregated from hourly data.
#' The original data spans from July 1, 2015, to June 30, 2024, and records total hourly electricity demand.
#' For this dataset, the hourly values have been aggregated by averaging across all 24 hours of each day,
#' resulting in a daily time series of average hourly electricity demand (measured in megawatts).
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Data.Date}{Date of the observation (class \code{Date}).}
#'   \item{Demand..MW.}{Daily average hourly electricity demand in megawatts (numeric).}
#' }
#'
#' @source \url{https://www.eia.gov/electricity/gridmonitor/dashboard/electric_overview/regional/REG-NY}
#' @references
#' U.S. Energy Information Administration (EIA). (2024). New York Independent System Operator electricity demand data.
#' Retrieved June 16, 2024, from \url{https://www.eia.gov/electricity/gridmonitor/dashboard/electric_overview/regional/REG-NY}.
#'
#' @examples
#' data(NYelectricity)
#' head(NYelectricity)
"NYelectricity"
