

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
When the BIC values are plotted, a number indicates the most probable number of ``alive'' clusters. The pairwise scatterplots (\code{what = "classification_pairs"}) are created using the \code{coordProj} function of the \code{mclust} package. The \code{what = "correlation"} is plotted using the \code{corrplot} package. Note that the \code{what = "classification_matplot"} plots the original data (before scaling and centering). On the other hand, the option \code{what = "classification_pairs"} plots the centered and scaled data. 
}
\author{
Panagiotis Papastamoulis
}
\references{
Luca Scrucca and Michael Fop and Thomas Brendan Murphy and Adrian E. Raftery (2017). mclust 5: clustering, classification and density estimation using Gaussian finite mixture models. The R Journal, 8(1): 205--233.

 Taiyun Wei and Viliam Simko (2017). R package "corrplot":
  Visualization of a Correlation Matrix (Version 0.84). Available from
  https://github.com/taiyun/corrplot

}
