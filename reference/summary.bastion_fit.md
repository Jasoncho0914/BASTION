# Component summaries of a BASTION fit

Returns posterior component means for fixed-period and calendar BASTION
fits.

## Usage

``` r
# S3 method for class 'bastion_fit'
summary(object, ...)
```

## Arguments

- object:

  A fitted object of class `bastion_fit`.

- ...:

  Additional arguments, currently unused.

## Value

The full T by M data frame, identical to `object$components`.
