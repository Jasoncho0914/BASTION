\name{plot.bastion_fit}
\alias{plot.bastion_fit}
\title{Plot a BASTION decomposition}
\usage{\method{plot}{bastion_fit}(x, ...)}
\arguments{\item{x}{A \code{bastion_fit} object.}\item{...}{Arguments passed to each base graphics plot.}}
\description{Draws observed data and posterior mean components with shaded credible bands in stacked panels on a common time axis. The interval level is taken from fitting metadata. Uses stored intervals and does not require retained draws.}
\value{The fitted object, invisibly.}
