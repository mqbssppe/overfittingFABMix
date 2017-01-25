\name{complete.log.likelihood}
\alias{complete.log.likelihood}
\title{
Complete log-likelihood function
}
\description{
Complete log-likelihood function
}
\usage{
complete.log.likelihood(x_data, w, mu, Lambda, SigmaINV, z)
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
  \item{Lambda}{
	Factor loadings
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


