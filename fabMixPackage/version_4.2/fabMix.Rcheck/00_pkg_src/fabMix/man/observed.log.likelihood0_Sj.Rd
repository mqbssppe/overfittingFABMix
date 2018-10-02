\name{observed.log.likelihood0_Sj}
\alias{observed.log.likelihood0_Sj}
\title{
Log-likelihood of the mixture model
}
\description{
Log-likelihood of the mixture model evaluated only at the alive components.
}
\usage{
observed.log.likelihood0_Sj(x_data, w, mu, Lambda, Sigma, z)
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
  \item{Lambda}{
	Factor loadings
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

