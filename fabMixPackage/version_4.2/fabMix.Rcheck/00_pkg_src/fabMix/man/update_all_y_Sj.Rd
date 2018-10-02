\name{update_all_y_Sj}
\alias{update_all_y_Sj}
\title{
Gibbs sampling for \eqn{y}
}
\description{
Gibbs sampling for \eqn{y}
}
\usage{
update_all_y_Sj(x_data, mu, SigmaINV, Lambda, z)
}
\arguments{
  \item{x_data}{
	Data
}
  \item{mu}{
	Marginal means
}
  \item{SigmaINV}{
	Precision matrix per component
}
  \item{Lambda}{
	Factor loadings
}
  \item{z}{
	Allocation vector
}
}
\value{
	A matrix with generated factors
}
\author{
Panagiotis Papastamoulis
}

