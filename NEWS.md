# BASTION 0.2.0

## New features and improvements

- Added calendar-aware decomposition with `fit_BASTION_calendar()`.
- Updated seasonality modeling with shared recurring patterns and locally adaptive deviations.
- Added `D_trend` to select first- or second-difference trend priors; the default is second differences.
- Restructured fitted objects for direct access to components, credible intervals, and recurring seasonal patterns.
- Added decomposition plots, component plots with credible intervals, and print and summary methods.
- Added simulation helpers for fixed-period and calendar-based examples.
- Added square-root QR fallback for Gaussian sampling.
- Updated the README, function documentation, and public reference pages.

## Bug fixes

- Fixed the outlier-model implementation to use the horseshoe+ prior.
- Removed the restriction fixing the first four outliers to zero and updated the observation-variance parameter count accordingly.
- Corrected the outlier column label when regression and outliers are both included.
- Corrected calendar plot axes to use the supplied timestamps' time zone.

# BASTION 0.1.0

- Initial BASTION R package.
