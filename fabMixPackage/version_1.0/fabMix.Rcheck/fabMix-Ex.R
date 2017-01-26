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
### Title: Overfitting Bayesian Mixtures of Factor Analyzers with an
###   Unknown Number of Components
### Aliases: fabMix-package
### Keywords: package

### ** Examples

# simulate a synthetic dataset along the lines of the paper:
n = 1000              # sample size
p = 40                # number of variables
q = 4                 # number of factors
K = 10                # number of clusters
sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
set.seed(10)
syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )

## Not run: 
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
##D dealWithLabelSwitching_same_sigma(x_data = syntheticDataset$data, 
##D         outputFolder = outputFolder, q = q, 
##D         compute_regularized_expression = TRUE, Km = Kmax)
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fabMix-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
