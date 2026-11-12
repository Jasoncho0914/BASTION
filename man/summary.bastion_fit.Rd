\name{summary.bastion_fit}
\alias{summary.bastion_fit}
\title{Component summaries of a BASTION fit}
\usage{
\method{summary}{bastion_fit}(object, ...)
}
\arguments{
\item{object}{A fitted object of class \code{bastion_fit}.}
\item{...}{Additional arguments, currently unused.}
}
\description{Returns posterior component means for fixed-period and calendar BASTION fits.}
\value{The full T by M data frame, identical to \code{object$components}.}
