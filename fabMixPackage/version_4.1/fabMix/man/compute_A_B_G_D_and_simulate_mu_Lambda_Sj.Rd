\name{compute_A_B_G_D_and_simulate_mu_Lambda_Sj}
\alias{compute_A_B_G_D_and_simulate_mu_Lambda_Sj}
\title{
	Computation and simulations
}
\description{
	This function simulates \eqn{\mu} and \eqn{\Lambda}.
}
\usage{
	compute_A_B_G_D_and_simulate_mu_Lambda_Sj(SigmaINV, 
		suff_statistics, OmegaINV, K, priorConst1, T_INV, v_r)
}
\arguments{
  \item{SigmaINV}{
	Precision matrix \eqn{\Sigma^{-1}} per component
}
  \item{suff_statistics}{
	Sufficient statistics
}
  \item{OmegaINV}{
	Prior parameter: \eqn{\Omega^{-1}}
}
  \item{K}{
	Number of overfitting mixture components
}
  \item{priorConst1}{
	Prior constant: \eqn{T^{-1}\xi}
}
  \item{T_INV}{
	Prior parameter: \eqn{T^{-1}\xi}
}
  \item{v_r}{
	This vector is used to set to zero the upper right \eqn{(q-1)\times(q-1)} diagonal block of factor loadings for identifiability purposes.
}
}
\value{
	A list containing \eqn{A}, \eqn{B}, \eqn{\Gamma}, \eqn{\Delta} and a draw from the conditional distributions of \eqn{\mu} and \eqn{\Lambda}.
}
\author{
Panagiotis Papastamoulis
}


