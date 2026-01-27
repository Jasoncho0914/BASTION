# Daily Average Electricity Demand

A dataset containing daily average electricity demand in New York State,
aggregated from hourly data. The original data spans from July 1, 2015,
to June 30, 2024, and records total hourly electricity demand. For this
dataset, the hourly values have been aggregated by averaging across all
24 hours of each day, resulting in a daily time series of average hourly
electricity demand (measured in megawatts).

## Usage

``` r
NYelectricity
```

## Format

A data frame with the following columns:

- Data.Date:

  Date of the observation (class `Date`).

- Demand..MW.:

  Daily average hourly electricity demand in megawatts (numeric).

## Source

<https://www.eia.gov/electricity/gridmonitor/dashboard/electric_overview/regional/REG-NY>

## References

U.S. Energy Information Administration (EIA). (2024). New York
Independent System Operator electricity demand data. Retrieved June 16,
2024, from
<https://www.eia.gov/electricity/gridmonitor/dashboard/electric_overview/regional/REG-NY>.

## Examples

``` r
data(NYelectricity)
head(NYelectricity)
#>    Data.Date Demand..MW.
#> 1 2015-07-01    20170.75
#> 2 2015-07-02    19227.46
#> 3 2015-07-03    17665.83
#> 4 2015-07-04    16489.88
#> 5 2015-07-05    17702.96
#> 6 2015-07-06    21223.83
```
