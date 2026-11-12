# Plot a BASTION decomposition

Draws observed data and posterior mean components with shaded credible
bands in stacked panels on a common time axis. The interval level is
taken from fitting metadata. Uses stored intervals and does not require
retained draws.

## Usage

``` r
# S3 method for class 'bastion_fit'
plot(x, ...)
```

## Arguments

- x:

  A `bastion_fit` object.

- ...:

  Arguments passed to each base graphics plot.

## Value

The fitted object, invisibly.
