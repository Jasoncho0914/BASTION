# Decomposition of time series data with Bayesian Adaptive Seasonality Trend decomposition Incorporating Outlier and Noise (BASTION)

decompose time series data using BASTION. The time series should be a
vector.

## Usage

``` r
fit_BASTION(
  y,
  Ks,
  X = NULL,
  Outlier = FALSE,
  cl = 0.95,
  sparse = FALSE,
  obsSV = "const",
  nchains = 2,
  nsave = 1000,
  nburn = 1000,
  nskip = 4,
  verbose = TRUE,
  save_samples = FALSE
)
```

## Arguments

- y:

  numeric vector of the `T x 1` vector of time series observations

- Ks:

  list of values containing the seasonal periods

- X:

  matrix for additional covariate for regression (default is NULL)

- Outlier:

  logical flag (default is FALSE) to model outlier

- cl:

  scalar between (0,1) for confidence leve (default is 0.95)

- sparse:

  logical flag (default is FALSE) to induce additional shrinkage to the
  trend estimate

- obsSV:

  options for modeling the error variance. It must be one of the
  following:

  - const: Constant error variance for all time points.

  - SV: Stochastic Volatility model.

- nchains:

  integer scalar for the number of chains for MCMC sampling (default is
  2)

- nsave:

  integer scalar (default = 1000); number of MCMC iterations to record

- nburn:

  integer scalar (default = 1000); number of MCMC iterations to discard
  (burn-in)

- nskip:

  integer scalar (default = 4); number of MCMC iterations to skip
  between saving iterations, i.e., save every (nskip + 1)th draw

- verbose:

  logical; report extra information on progress if true. Defaults to
  FALSE

- save_samples:

  logical; save and return posterior samples of each components if true.
  Default is FALSE

## Value

`fit_BASTION` returns an object class list.

## `summary`

A list providing summarized posterior estimates:

- `p_means`: A matrix containing the posterior means of the observed
  data (`y`) and each decomposed component: **Trend**, **Seasonality**,
  and **Outlier**, if applicable.

- `Trend_sum`: A matrix containing the posterior mean and 95% credible
  interval of the trend estimate.

- `Seasonal"k"_sum`: A matrix containing the posterior mean and 95%
  credible interval of the seasonal estimate, where `"k"` is determined
  by the input.

- `Outlier_sum`: A matrix containing the posterior mean and 95% credible
  interval of the outlier component.

## `samples`

(returned only when `save_samples` = TRUE) A list containing MCMC
samples of the relevant parameters:

- `beta_combined`: Posterior samples of **Trend + Seasonality**.

- `beta`: Posterior samples of each component excluding the remainder.

- `obs_sigma_t2`: Posterior samples of the variance of the observation
  equation.

- `evol_sigma_t2`: Posterior samples of the evolution error term.

- `remainder`: Posterior samples of the remainder term.

- `yhat`: Posterior samples of the signal + error term.
