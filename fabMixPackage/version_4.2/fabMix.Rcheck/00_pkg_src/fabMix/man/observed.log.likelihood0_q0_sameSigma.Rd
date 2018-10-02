\name{observed.log.likelihood0_q0_sameSigma}
\alias{observed.log.likelihood0_q0_sameSigma}
\title{
Log-likelihood of the mixture model for \eqn{q=0} and same variance of errors
}
\description{
Log-likelihood of the mixture model evaluated only at the alive components.
}
\usage{
observed.log.likelihood0_q0_sameSigma(x_data, w, mu, Sigma, z)
}
\arguments{
  \item{x_data}{
The observed data
}
  \item{w}{
	Vector of mixture weights
}
  \item{mu}{
	Vector of marginal means
}
  \item{Sigma}{
	Covariance matrix of the errors per cluster
}
  \item{z}{
	Allocation vector
}
}
\value{
	Log-likelihood value
}
\author{
Panagiotis Papastamoulis
}

