\name{update_z4_Sj}
\alias{update_z4_Sj}
\title{
	Collapsed Gibbs for \eqn{z}
}
\description{
	Collapsed Gibbs for \eqn{z}
}
\usage{
update_z4_Sj(w, mu, Lambda, SigmaINV, K, x_data)
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

