\name{overfittingMFA_UCC}
\alias{overfittingMFA_UCC}
\title{
	Basic MCMC sampler for CCC
}
\description{
	Gibbs sampling for fitting a UCC mixture model of factor analyzers.
}
\usage{
overfittingMFA_UCC(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, 
	g, h, alpha_prior, alpha_sigma, beta_sigma, 
	start_values, q, zStart, gibbs_z)
}
\arguments{
  \item{x_data}{
	normalized data
}
  \item{originalX}{
	observed raw data (only for plotting purpose)
}
  \item{outputDirectory}{
	Name of the output folder
}
  \item{Kmax}{
	Number of mixture components
}
  \item{m}{
	Number of iterations
}
  \item{thinning}{
	Thinning of chain
}
  \item{burn}{
	Burn-in period
}
  \item{g}{
	Prior parameter \eqn{g}. Default value: \eqn{g = 2}.
}
  \item{h}{
	Prior parameter \eqn{h}. Default value: \eqn{h = 1}.
}
  \item{alpha_prior}{
	Parameters of the Dirichlet prior distribution of mixture weights.
}
  \item{alpha_sigma}{
	Prior parameter \eqn{\alpha}. Default value: \eqn{\alpha = 2}.	
}
  \item{beta_sigma}{
	Prior parameter \eqn{\beta}. Default value: \eqn{\beta = 1}.
}
  \item{start_values}{
	Optional (not used)
}
  \item{q}{
	Number of factors.
}
  \item{zStart}{
	Optional (not used)
}
  \item{gibbs_z}{
	Optional
}
}
\value{
	List of files}
\author{
Panagiotis Papastamoulis
}


