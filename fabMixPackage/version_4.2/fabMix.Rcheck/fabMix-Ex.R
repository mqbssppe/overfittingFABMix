pkgname <- "fabMix"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "fabMix-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('fabMix')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("fabMix-package")
### * fabMix-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fabMix-package
### Title: Overfitting Parsimonious Bayesian Mixtures of Factor Analyzers
###   with an Unknown Number of Components
### Aliases: fabMix-package
### Keywords: package

### ** Examples


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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fabMix-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fabMix")
### * fabMix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fabMix
### Title: Main function
### Aliases: fabMix

### ** Examples

## Not run: 
##D # simulate a synthetic dataset along the lines of the paper:
##D n = 1000              # sample size
##D p = 40                # number of variables
##D q = 4                 # number of factors
##D K = 10                # number of clusters
##D sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
##D set.seed(10)
##D syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )
##D 
##D # define parameters
##D Kmax <- 20    # number of overfitted mixture components
##D nChains <- 8  # number of parallel chains
##D outputFolder <- "fabMixExample"
##D # Run algorithm
##D fabMix( rawData = syntheticDataset$data, 
##D 	nChains = nChains,
##D         outDir = outputFolder, Kmax = Kmax, 
##D 	mCycles = 1200, 
##D         burnCycles = 200, q = 3:5) 
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fabMix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fabMix_CxC")
### * fabMix_CxC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fabMix_CxC
### Title: Main function of the package for CUC, CCC models
### Aliases: fabMix_CxC

### ** Examples

## Not run: 
##D # simulate a synthetic dataset along the lines of the paper:
##D n = 1000              # sample size
##D p = 40                # number of variables
##D q = 4                 # number of factors
##D K = 10                # number of clusters
##D sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
##D set.seed(10)
##D syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )
##D 
##D # define parameters
##D Kmax <- 20    # number of overfitted mixture components
##D nChains <- 8  # number of parallel chains
##D dN <- 1
##D # Dirichlet prior of mixture weights per chain.
##D #   The target chain corresponds to the first entry.
##D dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax   
##D outputFolder <- "fabMixExample"
##D # Run algorithm
##D fabMix_CxC( dirPriorAlphas = dirPriorAlphas, 
##D         rawData = syntheticDataset$data, 
##D         outDir = outputFolder, Kmax = Kmax, mCycles = 1200, 
##D         burnCycles = 200, q = q) 
##D 
##D # Compute information criteria:
##D #getStuffForDIC(x_data = syntheticDataset$data, outputFolder = outputFolder, q = q)
##D 
##D # Deal with label switching:
##D #dealWithLabelSwitching(x_data = syntheticDataset$data, 
##D #        outputFolder = outputFolder, q = q, 
##D #        compute_regularized_expression = TRUE, Km = Kmax)
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fabMix_CxC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fabMix_CxU")
### * fabMix_CxU

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fabMix_CxU
### Title: Main function of the package for CCU, CUU models
### Aliases: fabMix_CxU

### ** Examples

## Not run: 
##D # simulate a synthetic dataset along the lines of the paper:
##D n = 1000              # sample size
##D p = 40                # number of variables
##D q = 4                 # number of factors
##D K = 10                # number of clusters
##D sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
##D set.seed(10)
##D syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )
##D 
##D # define parameters
##D Kmax <- 20    # number of overfitted mixture components
##D nChains <- 8  # number of parallel chains
##D dN <- 1
##D # Dirichlet prior of mixture weights per chain.
##D #   The target chain corresponds to the first entry.
##D dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax   
##D outputFolder <- "fabMixExample"
##D # Run algorithm
##D fabMix_CxU( dirPriorAlphas = dirPriorAlphas, 
##D         rawData = syntheticDataset$data, 
##D         outDir = outputFolder, Kmax = Kmax, mCycles = 1200, 
##D         burnCycles = 200, q = q) 
##D 
##D # Compute information criteria:
##D #getStuffForDIC(x_data = syntheticDataset$data, outputFolder = outputFolder, q = q)
##D 
##D # Deal with label switching:
##D #dealWithLabelSwitching(x_data = syntheticDataset$data, 
##D #        outputFolder = outputFolder, q = q, 
##D #        compute_regularized_expression = TRUE, Km = Kmax)
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fabMix_CxU", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fabMix_UxC")
### * fabMix_UxC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fabMix_UxC
### Title: Main function of the package for UUC, UCC models
### Aliases: fabMix_UxC

### ** Examples

## Not run: 
##D # simulate a synthetic dataset along the lines of the paper:
##D n = 1000              # sample size
##D p = 40                # number of variables
##D q = 4                 # number of factors
##D K = 10                # number of clusters
##D sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
##D set.seed(10)
##D syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )
##D 
##D # define parameters
##D Kmax <- 20    # number of overfitted mixture components
##D nChains <- 8  # number of parallel chains
##D dN <- 1
##D # Dirichlet prior of mixture weights per chain.
##D #   The target chain corresponds to the first entry.
##D dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax   
##D outputFolder <- "fabMixExample"
##D # Run algorithm
##D fabMix_UxC( dirPriorAlphas = dirPriorAlphas, 
##D         rawData = syntheticDataset$data, 
##D         outDir = outputFolder, Kmax = Kmax, mCycles = 1200, 
##D         burnCycles = 200, q = q) 
##D 
##D # Compute information criteria:
##D #getStuffForDIC(x_data = syntheticDataset$data, outputFolder = outputFolder, q = q)
##D 
##D # Deal with label switching:
##D #dealWithLabelSwitching(x_data = syntheticDataset$data, 
##D #        outputFolder = outputFolder, q = q, 
##D #        compute_regularized_expression = TRUE, Km = Kmax)
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fabMix_UxC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fabMix_UxU")
### * fabMix_UxU

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fabMix_UxU
### Title: Fit UxU mixtures
### Aliases: fabMix_UxU

### ** Examples

## Not run: 
##D # simulate a synthetic dataset along the lines of the paper:
##D n = 1000              # sample size
##D p = 40                # number of variables
##D q = 4                 # number of factors
##D K = 10                # number of clusters
##D sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
##D set.seed(10)
##D syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )
##D 
##D # define parameters
##D Kmax <- 20    # number of overfitted mixture components
##D nChains <- 8  # number of parallel chains
##D dN <- 1
##D # Dirichlet prior of mixture weights per chain.
##D #   The target chain corresponds to the first entry.
##D dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax   
##D outputFolder <- "fabMixExample"
##D # Run algorithm
##D fabMix_UxU( dirPriorAlphas = dirPriorAlphas, 
##D         rawData = syntheticDataset$data, 
##D         outDir = outputFolder, Kmax = Kmax, mCycles = 1200, 
##D         burnCycles = 200, q = q) 
##D 
##D # Compute information criteria:
##D getStuffForDIC(x_data = syntheticDataset$data, outputFolder = outputFolder, q = q)
##D 
##D # Deal with label switching:
##D dealWithLabelSwitching(x_data = syntheticDataset$data, 
##D         outputFolder = outputFolder, q = q, 
##D         compute_regularized_expression = TRUE, Km = Kmax)
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fabMix_UxU", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fabMix_missing_values")
### * fabMix_missing_values

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fabMix_missing_values
### Title: Main function for the case of missing values
### Aliases: fabMix_missing_values

### ** Examples

## Not run: 
##D # simulate a synthetic dataset along the lines of the paper:
##D n = 1000              # sample size
##D p = 40                # number of variables
##D q = 4                 # number of factors
##D K = 10                # number of clusters
##D sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
##D set.seed(10)
##D syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )
##D 
##D # define parameters
##D Kmax <- 20    # number of overfitted mixture components
##D nChains <- 8  # number of parallel chains
##D dN <- 1
##D # Dirichlet prior of mixture weights per chain.
##D #   The target chain corresponds to the first entry.
##D dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax   
##D outputFolder <- "fabMixExample"
##D # Run algorithm
##D fabMix( dirPriorAlphas = dirPriorAlphas, 
##D         rawData = syntheticDataset$data, 
##D         outDir = outputFolder, Kmax = Kmax, mCycles = 1200, 
##D         burnCycles = 200, q = q) 
##D 
##D # Compute information criteria:
##D getStuffForDIC(x_data = syntheticDataset$data, outputFolder = outputFolder, q = q)
##D 
##D # Deal with label switching:
##D dealWithLabelSwitching(x_data = syntheticDataset$data, 
##D         outputFolder = outputFolder, q = q, 
##D         compute_regularized_expression = TRUE, Km = Kmax)
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fabMix_missing_values", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
