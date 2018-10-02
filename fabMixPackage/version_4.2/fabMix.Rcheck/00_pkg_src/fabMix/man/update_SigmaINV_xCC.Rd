\name{update_SigmaINV_xCC}
\alias{update_SigmaINV_xCC}
\title{
Gibbs sampling for \eqn{\Sigma^{-1}} for xCC models
}
\description{
Gibbs sampling for \eqn{\Sigma^{-1}} for xCC models
}
\usage{
update_SigmaINV_xCC(x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma)
}
\arguments{
  \item{x_data}{
	Data
}
  \item{z}{
	Allocation vector
}
  \item{y}{
	Factors
}
  \item{Lambda}{
	Factor loadings
}
  \item{mu}{
	Marginal means
}
  \item{K}{
	Number of components
}
  \item{alpha_sigma}{
	Prior parameter
}
  \item{beta_sigma}{
	Prior parameter
}
}
\value{
	\eqn{\Sigma^{-1}}
}
\author{
Panagiotis Papastamoulis	
}

