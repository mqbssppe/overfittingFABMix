\name{fabMix}
\alias{fabMix}
\title{
Main function
}
\description{
This function runs parallel chains under a prior tempering scheme of the Dirichlet prior distribution of mixture weights. 
}
\usage{
fabMix(model, nChains, dirPriorAlphas, rawData, outDir, Kmax, mCycles, 
	burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize,  
	thinning, zStart, nIterPerCycle, gibbs_z, 
	warm_up_overfitting, warm_up, overfittingInitialization, 
	progressGraphs, gwar, rmDir)
}
\arguments{
  \item{model}{
Any subset of "UUU" "CUU" "UCU" "CCU" "UCC" "UUC" "CUC", "CCC" indicating the fitted models. By default, all models are fitted. 
}
  \item{nChains}{
The number of parallel heated chains. When `dirPriorAlphas` is supplied, `nChains` can be ignored.
}
  \item{dirPriorAlphas}{
vector of length \code{nChains} in the form of an increasing sequence of positive scalars. Each entry contains the (common) prior Dirichlet parameter for the corresponding chain. Default: \code{dirPriorAlphas = c(1, 1 + dN*(2:nChains - 1))/Kmax}, where \code{dN = 1}, for \code{nChains > 1}. Otherwise, \code{dirPriorAlphas = 1/Kmax}. 
}
  \item{rawData}{
The observed data as an \eqn{n\times p} matrix. Clustering is performed on the rows of the matrix.
}
  \item{outDir}{
Name of the output folder. An error is thrown if this directory already exists.
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
	Prior parameter \eqn{g}. Default value: \eqn{g = 0.5}.
}
  \item{h}{
	Prior parameter \eqn{h}. Default value: \eqn{h = 0.5}.
}
  \item{alpha_sigma}{
	Prior parameter \eqn{\alpha}. Default value: \eqn{\alpha = 0.5}.
}
  \item{beta_sigma}{
	Prior parameter \eqn{\beta}. Default value: \eqn{\beta = 0.5}.
}
  \item{q}{
	A vector containing the number of factors to be fitted.
}
  \item{normalize}{
	Should the observed data be normalized? Default value: TRUE. (Recommended) 
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
	Number of iterations for the overfitting initialization scheme. Default value: 500.
}
  \item{warm_up}{
	Number of iterations that will be used to initialize the models before starting proposing switchings. Default value: 5000.
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
  \item{rmDir}{
	Logical value indicating whether to delete the \code{outDir} directory. Default: TRUE.
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
\details{
Let \eqn{\boldsymbol{X}_i}; \eqn{i =1,\ldots,n} denote independent \eqn{p}-dimensional random vectors. Let \eqn{Y_i\in  R^q} with \eqn{q < p} denote the latent factor for observation \eqn{i = 1,\ldots,n}. In the typical factor analysis model, each observation is modelled as \eqn{\boldsymbol{X}_i = \boldsymbol{\mu} + \boldsymbol{\Lambda} \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i }, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\boldsymbol{\Sigma})}, where \eqn{\boldsymbol{\varepsilon}_i} and \eqn{Y_i} are assumed independent, \eqn{i =1,\ldots,n}. The \eqn{p\times q} matrix \eqn{\Lambda} consists of the factor loadings. Assume that there are \eqn{K} clusters and let \eqn{Z_i} denotes the latent allocation of observation \eqn{i} to one amongs the  \eqn{k} clusters, with prior probability \eqn{P(Z_i = k) = w_k}, \eqn{k = 1,\ldots,K}, independent for \eqn{i=1,\ldots,n}.  Following McNicholas et al (2008), the following parameterizations are used:

UUU model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda}_k \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\boldsymbol{\Sigma}_k)}

UCU model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda}_k \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\boldsymbol{\Sigma})}

UUC model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda}_k \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\sigma_k \boldsymbol{I}_p)}

UCC model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda}_k \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\sigma \boldsymbol{I}_p)}

CUU model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda} \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\boldsymbol{\Sigma}_k)}

CCU model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda} \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\boldsymbol{\Sigma})}

CUC model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda} \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\sigma_k \boldsymbol{I}_p)}

CCC model: \eqn{(\boldsymbol{X}_i|Z_i = k) = \boldsymbol{\mu}_k + \boldsymbol{\Lambda} \boldsymbol{Y}_i + \boldsymbol{\varepsilon}_i}, with \eqn{\boldsymbol{\varepsilon}_i \sim \mathcal N(0,\sigma \boldsymbol{I}_p)}

In all cases, \eqn{\boldsymbol{\varepsilon}_i} and \eqn{\boldsymbol{Y}_i} are assumed independent, \eqn{i =1,\ldots,n}. Note that \eqn{\boldsymbol{\Sigma}_k} and \eqn{\boldsymbol{\Sigma}} denote positive definite matrices, \eqn{\boldsymbol{I}_p} denotes the \eqn{p\times p} identity matrix and \eqn{\sigma_k}, \eqn{\sigma} denote positive scalars. 

}
\note{
It is recommended to always use: \code{normalize = TRUE} (default). Tuning of \code{dirPriorAlphas} may be necessary to achieve reasonable acceptance rates of chain swaps. Note that the output is reordered in order to deal with the label switching problem, according to the ECR algorithm applied by \code{dealWithLabelSwitching} function.
}
\seealso{
\code{\link{plot.fabMix.object}}
}
\references{
McNicholas, P.D. and Murphy, T.B. Stat Comput (2008) 18: 285. https://doi.org/10.1007/s11222-008-9056-0.

Papastamoulis, P. (2018). Overfitting Bayesian mixtures of factor analyzers with an unknown number of components. Computational Statistics and Data Analysis, 124: 220-234. DOI: 10.1016/j.csda.2018.03.007.
}
\author{
Panagiotis Papastamoulis
}
\examples{
\dontrun{
library('fabMix')

n = 90                # sample size
p = 50                # number of variables
q = 2                # number of factors
K = 3		     # number of clusters

sINV_diag = 1/((1:p))	 				# diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=FALSE,K.true = K, n = n, q = q, p = p, 
			sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
qRange <- 1:3	# range of values for the number of factors

Kmax <- 20		# number of components for the overfitted mixture model
nChains <- 8		# number of parallel heated chains

set.seed(3)
fm <- fabMix( rawData = syntheticDataset$data,
	nChains = nChains, outDir = "fabMixExample",
        Kmax = Kmax, mCycles = 500, burnCycles = 100, q = qRange) 


# Now print a run summary and produce some plots. 
print(fm)
plot(fm, what = "BIC")
plot(fm, what = "classification_pairs")
}


}


