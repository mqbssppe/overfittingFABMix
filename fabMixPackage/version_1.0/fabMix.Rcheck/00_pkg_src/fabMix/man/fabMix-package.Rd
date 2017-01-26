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

Papastamoulis P. and Iliopoulos G. (2010). An artificial allocations based solution to the label switching problem in Bayesian analysis of mixtures of distributions. Journal of Computational and Graphical Statistics, 19: 313-331.

van Havre, Z., White, N., Rousseau, J. and Mengersen, K. (2015). Overfitting Bayesian Mixture Models with an Unknown Number of Components. PLOS ONE, 10(7): 1-27.

Papastamoulis, P. (2016). \code{label.switching}: An R Package for Dealing with the Label Switching Problem in MCMC Outputs. Journal of Statistical Software, 69(1), 1-24.

Papastamoulis, P. (2017). Overfitting Bayesian mixtures of factor analyzers with an unknown number of components. arXiv:1701.04605 [stat.ME]
}
\keyword{ package }
\seealso{
\code{\link{fabMix}}, \code{\link{dealWithLabelSwitching_same_sigma}}, \code{\link{getStuffForDIC}}
}
\examples{
# simulate a synthetic dataset along the lines of the paper:
n = 1000              # sample size
p = 40                # number of variables
q = 4                 # number of factors
K = 10                # number of clusters
sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
set.seed(10)
syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )

\dontrun{

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

# Compute information criteria:
getStuffForDIC(x_data = syntheticDataset$data, outputFolder = outputFolder, q = q)

# Deal with label switching:
dealWithLabelSwitching_same_sigma(x_data = syntheticDataset$data, 
        outputFolder = outputFolder, q = q, 
        compute_regularized_expression = TRUE, Km = Kmax)
}


}
