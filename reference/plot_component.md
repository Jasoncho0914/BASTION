# Plot a BASTION component with credible intervals

Uses stored posterior means and credible intervals. Recurring patterns
use their native seasonal positions; other panels use observation times.
Does not require retained draws.

## Usage

``` r
plot_component(fit, component,
  part = c("all", "total", "recurring", "deviation"),
  col = "black", band_col = grDevices::adjustcolor("#0072B2", .50), ...)
```

## Arguments

- fit:

  A `bastion_fit` object.

- component:

  Exact column name in `fit$components`.

- part:

  Seasonal view. The default shows total, recurring, and deviation in
  separate panels.

- col:

  Posterior mean line color.

- band_col:

  Credible band fill color.

- ...:

  Arguments passed to each base graphics plot.

## Value

The fitted object, invisibly.
