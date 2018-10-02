\name{fabMix-package}
\alias{fabMix-package}
\docType{package}
\title{
\packageTitle{fabMix}
}
\description{
\packageDescription{fabMix}
}
\author{
\packageAuthor{fabMix}

Maintainer: \packageMaintainer{fabMix}
}
\references{
Fokoue, E. and Titterington, D.M. (2003). Mixtures of Factor Analysers: Bayesian Estimation and Inference by Stochastic Simulation. Machine Learing, 50(1): 73-94.

McNicholas, P.D. and Murphy, T.B. Stat Comput (2008) 18: 285. https://doi.org/10.1007/s11222-008-9056-0.

Papastamoulis P. and Iliopoulos G. (2010). An artificial allocations based solution to the label switching problem in Bayesian analysis of mixtures of distributions. Journal of Computational and Graphical Statistics, 19: 313-331.

van Havre, Z., White, N., Rousseau, J. and Mengersen, K. (2015). Overfitting Bayesian Mixture Models with an Unknown Number of Components. PLOS ONE, 10(7): 1-27.

Papastamoulis, P. (2016). \code{label.switching}: An R Package for Dealing with the Label Switching Problem in MCMC Outputs. Journal of Statistical Software, 69(1), 1-24.

Papastamoulis, P. (2018). Overfitting Bayesian mixtures of factor analyzers with an unknown number of components. Computational Statistics and Data Analysis, 124: 220-234. DOI: 10.1016/j.csda.2018.03.007.
}
\keyword{ package }
\seealso{
\code{\link{fabMix}}, \code{\link{plot.fabMix.object}}
}
\examples{

library('fabMix')

n = 10               # sample size
p = 8                # number of variables
q = 2                # number of factors
K = 2		     # number of clusters

sINV_diag = 1/((1:p))	 				# diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
			sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
qRange <- 1:2	# range of values for the number of factors

Kmax <- 5		# number of components for the overfitted mixture model
nChains <- 2		# number of parallel heated chains

# Run `fabMix` for a _small_ number of iterations for the 
#	`UUU` (maximal model) and `CCC` (minimal model) parameterizations,
# 	using the default prior parallel heating parameters `dirPriorAlphas`.
#	NOTE: `dirPriorAlphas` may require some tuning in general.

set.seed(1)
fm <- fabMix( model = c("UUU", "CCC"), nChains = 2, 
	rawData = syntheticDataset$data, outDir = "toyExample",
        Kmax = Kmax, mCycles = 4, burnCycles = 1, q = qRange,
        g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
        warm_up_overfitting = 5, warm_up = 25) 

# WARNING: the following parameters: 
#  Kmax, nChains, mCycles, burnCycles, warm_up_overfitting, warm_up 
#	 should take (much) _larger_ values. E.g. a typical implementation consists of:
#        Kmax = 20, nChains = 8, mCycles = 1100, burnCycles = 100, 
#        warm_up_overfitting = 500, warm_up = 5000. 

# Now print a run summary and produce some plots. 
print(fm)
plot(fm, what = "BIC")
plot(fm, what = "classification_pairs")

}
