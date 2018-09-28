\name{update_z4}
\alias{update_z4}
\title{
	Collapsed Gibbs for \eqn{z}
}
\description{
	Collapsed Gibbs for \eqn{z}
}
\usage{
update_z4(w, mu, Lambda, SigmaINV, K, x_data)
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
  \item{SigmaINV}{
	Precision matrix
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

