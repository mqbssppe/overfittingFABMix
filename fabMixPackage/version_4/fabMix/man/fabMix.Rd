\name{fabMix}
\alias{fabMix}
\title{
Main function
}
\description{
This function runs parallel chains under a prior tempering scheme of the Dirichlet prior distribution of mixture weights. 
}
\usage{
fabMix(model, dirPriorAlphas, rawData, outDir, Kmax, mCycles, 
	burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize,  
	thinning, zStart, nIterPerCycle, gibbs_z, 
	warm_up_overfitting, warm_up, overfittingInitialization, 
	progressGraphs, gwar)
}
\arguments{
  \item{model}{
Any subset of "UUU" "CUU" "UCU" "CCU" "UCC" "UUC" "CUC", "CCC" indicating the fitted models.
}
  \item{dirPriorAlphas}{
The prior Dirichlet parameters for each chain.
}
  \item{rawData}{
The observed data as an \eqn{n\times p} matrix. Clustering is performed on the rows of the matrix.
}
  \item{outDir}{
Name of the output folder.
}
  \item{Kmax}{
Number of components in the overfitted mixture. Default: 20.
}
  \item{mCycles}{
Number of MCMC cycles. Each cycle consists of \code{nIterPerCycle} MCMC iterations. At the end of each cycle a swap of the state of two randomly chosen adjacent chains is attempted.
}
  \item{burnCycles}{
Number of cycles that will be discarded as burn-in period.
}
  \item{g}{
	Prior parameter \eqn{g}. Default value: \eqn{g = 2}.
}
  \item{h}{
	Prior parameter \eqn{h}. Default value: \eqn{h = 1}.
}
  \item{alpha_sigma}{
	Prior parameter \eqn{\alpha}. Default value: \eqn{\alpha = 2}.
}
  \item{beta_sigma}{
	Prior parameter \eqn{\beta}. Default value: \eqn{\beta = 1}.
}
  \item{q}{
	A vector containing the number of factors to be fitted.
}
  \item{normalize}{
	Should the observed data be normalized? Default value: TRUE. 
}
  \item{thinning}{
	Optional integer denoting the thinning of the keeped MCMC cycles.
}
  \item{zStart}{
	Optional starting value for the allocation vector.
}
  \item{nIterPerCycle}{
	Number of iteration per MCMC cycle. Default value: 10.
}
  \item{gibbs_z}{
	Select the gibbs sampling scheme for updating latent allocations of mixture model. Default value: 1.
}
  \item{warm_up_overfitting}{
	Number of iterations for the overfitting initialization scheme. Default value: 100.
}
  \item{warm_up}{
	Number of iterations that will be used to initialize the models before starting proposing switchings. Default value: 500.
}
  \item{overfittingInitialization}{
	Logical value indicating whether the chains are initialized via the overfitting initialization scheme. Default: TRUE.
}
  \item{progressGraphs}{
	Logical value indicating whether to plot successive states of the chains while the sampler runs. Default: FALSE.
}
  \item{gwar}{
	Initialization parameter. Default: 0.05.
}
}
\value{
	An object of class \code{fabMix.object}, that is, a list consisting of the following entries:
  \item{bic }{Bayesian Information Criterion per model and number of factors.}
  \item{class }{The estimated single best clustering of the observations according to the selected model.}
  \item{n_Clusters_per_model }{The most probable number of clusters (number of non-empty components of the overfitted mixture) per model and number of factors.}
  \item{posterior_probability }{The posterior probability of the estimated allocations according to the selected model.}
  \item{covariance_matrix }{The estimated posterior mean of the covariance matrix per cluster according to the selected model.}
  \item{mu }{The estimated posterior mean of the mean per cluster according to the selected model.}
  \item{weights }{The estimated posterior mean of the mixing proportions according to the selected model.}
  \item{selected_model }{Data frame containing the parameterization, number of clusters and factors of the selected model.}
  \item{mcmc }{A list containing the MCMC draws for the parameters of the selected models. All component-specific parameters have been reordered according to the ECR algorithm in order to undo the label switching problem.}
  \item{data }{The observed data.}
}
\note{
It is recommended to always use: \code{normalize = TRUE} (default). Tuning of \code{dirPriorAlphas} may be necessary to achieve reasonable acceptance rates of chain swaps. Note that the output is reordered in order to deal with the label switching problem, according to the ECR algorithm applied by \code{dealWithLabelSwitching} function.
}
\seealso{
\code{\link{dealWithLabelSwitching}}
}
\references{
Papastamoulis, P. (2017). Overfitting Bayesian mixtures of factor analyzers with an unknown number of components. arXiv:1701.04605 [stat.ME]
}
\author{
Panagiotis Papastamoulis
}
\examples{
\dontrun{
# simulate a synthetic dataset along the lines of the paper:
n = 1000              # sample size
p = 40                # number of variables
q = 4                 # number of factors
K = 10                # number of clusters
sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
set.seed(10)
syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )

# define parameters
Kmax <- 20    # number of overfitted mixture components
nChains <- 8  # number of parallel chains
dN <- 1
# Dirichlet prior of mixture weights per chain.
#   The target chain corresponds to the first entry.
dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax   
outputFolder <- "fabMixExample"
# Run algorithm
fabMix( dirPriorAlphas = dirPriorAlphas, 
        rawData = syntheticDataset$data, 
        outDir = outputFolder, Kmax = Kmax, mCycles = 1200, 
        burnCycles = 200, q = q) 
}


}


