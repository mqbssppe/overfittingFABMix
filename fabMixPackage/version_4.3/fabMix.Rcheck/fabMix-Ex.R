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

#################################################################
# (a) using 2 cores in parallel, each one running 2 heated chains.
#################################################################
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
fm <- fabMix( model = c("UUU", "CCC"), nChains = nChains, 
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

#################################################################
# (b) using 12 cores_____________________________________________
#_______4 models with 3 heated chains running in parallel________
#_______considering all 8 model parameterizations________________
#################################################################
## Not run: 
##D library('fabMix')
##D set.seed(99)
##D n = 100                # sample size
##D p = 30                # number of variables
##D q = 2                # number of factors
##D K = 5		     # number of clusters
##D sINV_diag = rep(1/100,p) 	# diagonal of inverse variance of errors
##D syntheticDataset <- simData(sameLambda=FALSE,K.true = K, n = n, q = q, p = p, 
##D 			sINV_values = sINV_diag)
##D colnames(syntheticDataset$data) <- paste0("x_",1:p)
##D qRange <- 1:3	# range of values for the number of factors
##D Kmax <- 20	# number of components for the overfitted mixture model
##D nChains <- 3	# number of parallel heated chains
##D 
##D fm <- fabMix( parallelModels = 4, 
##D 	nChains = nChains, 
##D 	model = c("UUU","CUU","UCU","CCU","UCC","UUC","CUC","CCC"), 
##D 	rawData = syntheticDataset$data, outDir = "toyExample_b",
##D         Kmax = Kmax, mCycles = 600, burnCycles = 100, q = qRange,
##D         g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
##D         warm_up_overfitting = 500, warm_up = 5000) 
##D print(fm)
##D plot(fm, what = "BIC")
##D plot(fm, what = "classification_pairs")
##D 
## End(Not run)




cleanEx()
nameEx("fabMix")
### * fabMix

flush(stderr()); flush(stdout())

### Name: fabMix
### Title: Main function
### Aliases: fabMix

### ** Examples

#################################################################
# (b) using 12 cores_____________________________________________
#_______4 models with 3 heated chains running in parallel________
#_______considering all 8 model parameterizations________________
#################################################################
## Not run: 
##D library('fabMix')
##D set.seed(99)
##D n = 100                # sample size
##D p = 30                # number of variables
##D q = 2                # number of factors
##D K = 5		     # number of clusters
##D sINV_diag = rep(1/100,p) 	# diagonal of inverse variance of errors
##D syntheticDataset <- simData(sameLambda=FALSE,K.true = K, n = n, q = q, p = p, 
##D 			sINV_values = sINV_diag)
##D colnames(syntheticDataset$data) <- paste0("x_",1:p)
##D qRange <- 1:3	# range of values for the number of factors
##D Kmax <- 20	# number of components for the overfitted mixture model
##D nChains <- 3	# number of parallel heated chains
##D 
##D fm <- fabMix( parallelModels = 4, 
##D 	nChains = nChains, 
##D 	model = c("UUU","CUU","UCU","CCU","UCC","UUC","CUC","CCC"), 
##D 	rawData = syntheticDataset$data, outDir = "toyExample_b",
##D         Kmax = Kmax, mCycles = 600, burnCycles = 100, q = qRange,
##D         g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
##D         warm_up_overfitting = 500, warm_up = 5000) 
##D print(fm)
##D plot(fm, what = "BIC")
##D plot(fm, what = "classification_pairs")
##D 
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
