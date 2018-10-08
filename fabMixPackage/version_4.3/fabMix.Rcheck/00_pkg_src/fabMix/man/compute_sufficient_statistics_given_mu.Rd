\name{compute_sufficient_statistics_given_mu}
\alias{compute_sufficient_statistics_given_mu}
\title{
Compute sufficient statistics given mu
}
\description{
Compute sufficient statistics given \eqn{y} and \eqn{z}.
}
\usage{
compute_sufficient_statistics_given_mu(y, z, K, x_data,mu)
}
\arguments{
  \item{y}{
	Matrix of factors
}
  \item{z}{
	Allocation vector
}
  \item{K}{
	Number of components
}
  \item{x_data}{
	Data
}
  \item{mu}{
	Means per component
}
}
\value{
	A list with six entries of sufficient statistics.
}
\author{
Panagiotis Papastamoulis
}

