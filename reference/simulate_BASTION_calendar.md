# Simulate data with known calendar seasonality

Generate observations, timestamps, and known components for calendar
fitting examples.

## Usage

``` r
simulate_BASTION_calendar(scenario = c("daily", "biweekly", "hourly"), seed = NULL)
```

## Arguments

- scenario:

  One of "daily", "biweekly", or "hourly"; default "daily".

- seed:

  Optional nonnegative integer random seed.

## Details

Each scenario generates y as the sum of its known components and
independent zero-mean Gaussian noise. Timestamps use UTC.

`"daily"` generates daily observations from 2019 through 2021, including
February 29, 2020. With t starting at zero, the trend is 50 + 0.005\*t.
Using midpoint calendar phases w within the week and a within the year,
weekly seasonality is 3\*cos(2\*pi\*w), and yearly seasonality is
8\*sin(2\*pi\*a) + 2\*cos(4\*pi\*a). Noise SD is 1.5.

`"biweekly"` generates 560 daily observations starting January 6, 2020.
With t starting at zero, the trend is 100 + 0.01\*t +
3\*sin(2\*pi\*t/560), and seasonality is 6\*sin(2\*pi\*t/14) +
2\*cos(4\*pi\*t/14). Noise SD is 1.5.

`"hourly"` generates 2016 hourly observations starting January 1, 2024
(12 full weeks). With t starting at zero, the trend is 100 + 0.002\*t.
Using midpoint phase h within each day, daily seasonality is
10\*sin(2\*pi\*h) + 2\*cos(4\*pi\*h). Monday-to-Sunday effects are c(2,
3, 3, 2, 1, -5, -6), constant within each day. Noise SD is 2.

## Value

A list with numeric observations `y`, POSIXct timestamps `dates`, a data
frame `truth` containing the trend and seasonal components, the Gaussian
`noise` vector, and its generating `noise_sd`.

## See also

[`fit_BASTION_calendar`](https://jasoncho0914.github.io/BASTION/reference/fit_BASTION_calendar.md)
