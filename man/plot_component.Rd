\name{plot_component}
\alias{plot_component}
\title{Plot a BASTION component with credible intervals}
\usage{plot_component(fit, component,
  part = c("all", "total", "recurring", "deviation"),
  col = "black", band_col = grDevices::adjustcolor("#0072B2", .50), ...)}
\arguments{
\item{fit}{A \code{bastion_fit} object.}
\item{component}{Exact column name in \code{fit$components}.}
\item{part}{Seasonal view. The default shows total, recurring, and deviation in separate panels.}
\item{col}{Posterior mean line color.}
\item{band_col}{Credible band fill color.}
\item{...}{Arguments passed to each base graphics plot.}}
\description{Uses stored posterior means and credible intervals. Recurring patterns use their native seasonal positions; other panels use observation times. Does not require retained draws.}
\value{The fitted object, invisibly.}
