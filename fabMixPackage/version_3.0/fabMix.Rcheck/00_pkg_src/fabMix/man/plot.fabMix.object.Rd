

\name{plot.fabMix.object}
\alias{plot.fabMix.object}
\title{
Plot function
}
\description{
This function plots \code{fabMix} function.
}
\usage{
\method{plot}{fabMix.object}(x, what, variableSubset, ...)
}
\arguments{
  \item{x}{
An object of class \code{fabMix.object}, which is returned by the \code{fabMix} function.
}
  \item{what}{
One of the "BIC", "classification_matplot", "classification_pairs", "correlation", "regularized_expression". The plot will display the BIC values per model and number of factors (along with the most probable number of clusters as text), a matplot per cluster for the selected model, scatterplots pairs, the estimated correlation matrix per cluster, and the estimated regularized expression of each variable to the factor space for the selected model, respectively.
}
  \item{variableSubset}{
An optional subset of the variables. By default, all variables are selected.
}
  \item{...}{
        ignored.
}
}
\details{
The function plots the BIC values obtained for each of the fitted models, visualizes the clusters.
}
\author{
Panagiotis Papastamoulis
}

