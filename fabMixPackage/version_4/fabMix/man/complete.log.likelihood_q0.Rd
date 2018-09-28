\name{complete.log.likelihood_q0}
\alias{complete.log.likelihood_q0}
\title{
Complete log-likelihood function for \eqn{q=0}
}
\description{
Complete log-likelihood function
}
\usage{
complete.log.likelihood_q0(x_data, w, mu, SigmaINV, z)
}

\arguments{
  \item{x_data}{
	Data
}
  \item{w}{
	Mixture weights
}
  \item{mu}{
	Marginal means
}
  \item{SigmaINV}{
	Precision matrix (inverse covariance)
}
  \item{z}{
	Allocation vector of the data to the mixture components
}
}
\value{
	complete log-likelihood value
}
\author{
Panagiotis Papastamoulis
}


