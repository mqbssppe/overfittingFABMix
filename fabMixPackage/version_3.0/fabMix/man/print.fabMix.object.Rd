

\name{print.fabMix.object}
\alias{print.fabMix.object}
\title{
Print function
}
\description{
This function prints a summary of objects returned by the \code{fabMix} function.
}
\usage{
\method{print}{fabMix.object}(x, printSubset, ...)
}
\arguments{
  \item{x}{
An object of class \code{fabMix.object}, which is returned by the \code{fabMix} function.
}
  \item{printSubset}{
Logical indicating whether to print the header or the whole matrix of estimates. Default value: TRUE.
}
  \item{...}{
        ignored.
}
}
\details{
The function prints the estimated distribution of the number of clusters, the estimated number of observations assigned to each cluster after post-processing the output with three label switching algorithms, as well as the header of the posterior mean estimates of \eqn{\theta_{kj}} (probability of success for cluster \eqn{k} and feature \eqn{j}) (conditionally on the most probable number of clusters).
}
\author{
Panagiotis Papastamoulis
}

