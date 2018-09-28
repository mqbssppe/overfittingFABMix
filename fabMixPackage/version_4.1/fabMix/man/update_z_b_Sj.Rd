\name{update_z_b_Sj}
\alias{update_z_b_Sj}
\title{
Gibbs sampling for \eqn{z}
}
\description{
Gibbs sampling for \eqn{z}
}
\usage{
update_z_b_Sj(w, mu, Lambda, y, SigmaINV, K, x_data)
}
\arguments{
  \item{w}{
	Mixture weights
}
  \item{mu}{
	Marginal means
}
  \item{Lambda}{
	Factor loadings
}
  \item{y}{
	Matrix of factors
}
  \item{SigmaINV}{
	Precision matrix per component
}
  \item{K}{
	Number of components
}
  \item{x_data}{
	Data
}
}
\value{
	Allocation vector
}
\author{
Panagiotis Papastamoulis
}

