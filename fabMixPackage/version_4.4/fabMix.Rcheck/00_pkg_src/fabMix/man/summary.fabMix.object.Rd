

\name{summary.fabMix.object}
\alias{summary.fabMix.object}
\title{
Summary method
}
\description{
S3 method for printing a summary of a \code{fabMix.object}.
}
\usage{
\method{summary}{fabMix.object}(object, ...)
}
\arguments{
  \item{object}{
An object of class \code{fabMix.object}, which is returned by the \code{fabMix} function.
}
  \item{...}{Ignored.}
}
\details{
The function prints a summary of the estimated posterior means for the parameters of the selected model for a \code{fabMix.object}. In particular, the method prints the ergodic means of the mixing proportions, marginal means and covariance matrix per component.  
}
\note{
More detailed summaries for the MCMC output can be computed by applying the \code{summary} function of the \code{coda} package to the \code{mcmc} object \code{object$mcmc}.
}
\author{
Panagiotis Papastamoulis
}

