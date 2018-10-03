pkgname <- "fabMix"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('fabMix')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("fabMix-package")
### * fabMix-package

flush(stderr()); flush(stdout())

### Name: fabMix-package
### Title: Overfitting Bayesian Mixtures of Factor Analyzers with
###   Parsimonious Covariance and Unknown Number of Components
### Aliases: fabMix-package
### Keywords: package

### ** Examples

# TOY EXAMPLE (very small numbers...)
library('fabMix')

n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2		     # number of clusters

sINV_diag = 1/((1:p))	 				# diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
			sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
qRange <- 1:2	# range of values for the number of factors

Kmax <- 20		# number of components for the overfitted mixture model
nChains <- 2		# number of parallel heated chains

# Run `fabMix` for a _small_ number of iterations for the 
#	`UUU` (maximal model) and `CCC` (minimal model) parameterizations,
# 	using the default prior parallel heating parameters `dirPriorAlphas`.
#	NOTE: `dirPriorAlphas` may require some tuning in general.

set.seed(3)
fm <- fabMix( model = c("UUU", "CCC"), nChains = 2, 
	rawData = syntheticDataset$data, outDir = "toyExample",
        Kmax = Kmax, mCycles = 4, burnCycles = 1, q = qRange,
        g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
        warm_up_overfitting = 2, warm_up = 3) 

# WARNING: the following parameters: 
#  nChains, mCycles, burnCycles, warm_up_overfitting, warm_up 
#	 should take (much) _larger_ values. E.g. a typical implementation consists of:
#        nChains = 8, mCycles = 1100, burnCycles = 100, 
#        warm_up_overfitting = 500, warm_up = 5000. 

# Now print a run summary and produce some plots. 
print(fm)
plot(fm, what = "BIC")
plot(fm, what = "classification_pairs")




cleanEx()
nameEx("fabMix")
### * fabMix

flush(stderr()); flush(stdout())

### Name: fabMix
### Title: Main function
### Aliases: fabMix

### ** Examples

## Not run: 
##D library('fabMix')
##D 
##D n = 90                # sample size
##D p = 50                # number of variables
##D q = 2                # number of factors
##D K = 3		     # number of clusters
##D 
##D sINV_diag = 1/((1:p))	 				# diagonal of inverse variance of errors
##D set.seed(100)
##D syntheticDataset <- simData(sameLambda=FALSE,K.true = K, n = n, q = q, p = p, 
##D 			sINV_values = sINV_diag)
##D colnames(syntheticDataset$data) <- paste0("x_",1:p)
##D qRange <- 1:3	# range of values for the number of factors
##D 
##D Kmax <- 20		# number of components for the overfitted mixture model
##D nChains <- 8		# number of parallel heated chains
##D 
##D set.seed(3)
##D fm <- fabMix( rawData = syntheticDataset$data,
##D 	nChains = nChains, outDir = "fabMixExample",
##D         Kmax = Kmax, mCycles = 500, burnCycles = 100, q = qRange) 
##D 
##D 
##D # Now print a run summary and produce some plots. 
##D print(fm)
##D plot(fm, what = "BIC")
##D plot(fm, what = "classification_pairs")
## End(Not run)





### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
