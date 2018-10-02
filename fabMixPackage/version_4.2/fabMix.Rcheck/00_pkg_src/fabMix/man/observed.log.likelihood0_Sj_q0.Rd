\name{observed.log.likelihood0_Sj_q0}
\alias{observed.log.likelihood0_Sj_q0}
\title{
Log-likelihood of the mixture model for \eqn{q=0}
}
\description{
Log-likelihood of the mixture model evaluated only at the alive components.
}
\usage{
observed.log.likelihood0_Sj_q0(x_data, w, mu, Sigma, z)
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

