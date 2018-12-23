pkgname <- "fabMix"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('fabMix')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("complete.log.likelihood")
### * complete.log.likelihood

flush(stderr()); flush(stdout())

### Name: complete.log.likelihood
### Title: Complete log-likelihood function for xCx models.
### Aliases: complete.log.likelihood

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- waveDataset1500[ 1:20, -1] # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( runif(K * p), dim = c(K,p) )
	Lambda <- array( runif(K*p*q), dim = c(K,p,q) )
	SigmaINV <- array(1, dim = c(p,p))
	# compute the complete.log.likelihood
	complete.log.likelihood(x_data = x_data, w = w, mu = mu, 
		Lambda = Lambda, SigmaINV = SigmaINV, z = z)



cleanEx()
nameEx("complete.log.likelihood_Sj")
### * complete.log.likelihood_Sj

flush(stderr()); flush(stdout())

### Name: complete.log.likelihood_Sj
### Title: Complete log-likelihood function for xUx models.
### Aliases: complete.log.likelihood_Sj

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- waveDataset1500[ 1:20, -1] # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( runif(K * p), dim = c(K,p) )
	Lambda <- array( runif(K*p*q), dim = c(K,p,q) )
	SigmaINV <- array( c(0.5, 0.75, 1), dim = c(K,p,p))
	# compute the complete.log.likelihood
	complete.log.likelihood_Sj(x_data = x_data, w = w, mu = mu, 
		Lambda = Lambda, SigmaINV = SigmaINV, z = z)



cleanEx()
nameEx("complete.log.likelihood_q0")
### * complete.log.likelihood_q0

flush(stderr()); flush(stdout())

### Name: complete.log.likelihood_q0
### Title: Complete log-likelihood function for xUx models and q=0
### Aliases: complete.log.likelihood_q0

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- scale(waveDataset1500[ 1:20, -1]) # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( -0.1 + 0.2*runif(K * p), dim = c(K,p) )
	SigmaINV <- array( 1, dim = c(K,p,p))
	# compute the complete.log.likelihood ( -inf )
	complete.log.likelihood_q0(x_data = x_data, w = w, mu = mu, 
		SigmaINV = SigmaINV, z = z)



cleanEx()
nameEx("complete.log.likelihood_q0_sameSigma")
### * complete.log.likelihood_q0_sameSigma

flush(stderr()); flush(stdout())

### Name: complete.log.likelihood_q0_sameSigma
### Title: Complete log-likelihood function for xCx models and q=0
### Aliases: complete.log.likelihood_q0_sameSigma

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- scale(waveDataset1500[ 1:20, -1]) # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( -0.1 + 0.2*runif(K * p), dim = c(K,p) )
	SigmaINV <- array( 1, dim = c(p,p))
	# compute the complete.log.likelihood ( -inf )
	complete.log.likelihood_q0_sameSigma(x_data = x_data, w = w, mu = mu, 
		SigmaINV = SigmaINV, z = z)



cleanEx()
nameEx("compute_A_B_G_D_and_simulate_mu_Lambda")
### * compute_A_B_G_D_and_simulate_mu_Lambda

flush(stderr()); flush(stdout())

### Name: compute_A_B_G_D_and_simulate_mu_Lambda
### Title: Computation and simulations
### Aliases: compute_A_B_G_D_and_simulate_mu_Lambda

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- scale(as.matrix(waveDataset1500[ 1:20, -1])) # data
	z <-  waveDataset1500[ 1:20, 1] # class
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	q <- 2
	K <- length(table(z))           # 3 classes

	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colMeans(x_data)
	priorConst1 <- array(diag(T_INV)*ksi, dim =c(p,1))
	# give some arbitrary values to the parameters:
	set.seed(1)
	mu <- array( runif(K * p), dim = c(K,p) )
	y <- array(rnorm(n = q*n), dim = c(n,q))
	SigmaINV <- array(data = 0, dim = c(p,p) )
	diag(SigmaINV) <- 0.5 + 0.5*runif(p)
	OmegaINV <- diag(q)
	# compute sufficient stats 
	suf_stat <- compute_sufficient_statistics(y = y, 
	 z = z, K = K, x_data = x_data)

	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}
	# now simulate mu and Lambda
	f2 <- compute_A_B_G_D_and_simulate_mu_Lambda(SigmaINV = SigmaINV, 
                suff_statistics = suf_stat, OmegaINV = OmegaINV, 
                K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
	# f2$mu contains the simulated means
	# f2$Lambdas contains the simulated factor loadings




cleanEx()
nameEx("compute_A_B_G_D_and_simulate_mu_Lambda_CCU")
### * compute_A_B_G_D_and_simulate_mu_Lambda_CCU

flush(stderr()); flush(stdout())

### Name: compute_A_B_G_D_and_simulate_mu_Lambda_CCU
### Title: Computation and simulations for CCU
### Aliases: compute_A_B_G_D_and_simulate_mu_Lambda_CCU

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- scale(as.matrix(waveDataset1500[ 1:20, -1])) # data
	z <-  waveDataset1500[ 1:20, 1] # class
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	q <- 2
	K <- length(table(z))           # 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	mu <- array( runif(K * p), dim = c(K,p) )
	y <- array(rnorm(n = q*n), dim = c(n,q))
	SigmaINV <- array(data = 0, dim = c(p,p) )
	diag(SigmaINV) = 0.5 + 0.5*runif(p)
	OmegaINV <- diag(q)
	# compute sufficient stats 
	suf_stat <- compute_sufficient_statistics_given_mu(y = y, 
	 z = z, K = K, x_data = x_data, mu = mu)

	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colMeans(x_data)
	priorConst1 <- array(diag(T_INV)*ksi, dim =c(p,1))
	# now simulate mu and Lambda
	f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_CCU(SigmaINV = SigmaINV, 
                suff_statistics = suf_stat, OmegaINV = OmegaINV, 
                K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
	# f2$mu contains the simulated means
	# f2$Lambdas contains the simulated factor loadings




cleanEx()
nameEx("compute_A_B_G_D_and_simulate_mu_Lambda_CUU")
### * compute_A_B_G_D_and_simulate_mu_Lambda_CUU

flush(stderr()); flush(stdout())

### Name: compute_A_B_G_D_and_simulate_mu_Lambda_CUU
### Title: Computation and simulations for CUU
### Aliases: compute_A_B_G_D_and_simulate_mu_Lambda_CUU

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- scale(as.matrix(waveDataset1500[ 1:20, -1])) # data
	z <-  waveDataset1500[ 1:20, 1] # class
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	q <- 2
	K <- length(table(z))           # 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	mu <- array( runif(K * p), dim = c(K,p) )
	y <- array(rnorm(n = q*n), dim = c(n,q))
	SigmaINV <- array(data = 0, dim = c(K,p,p) )
	for(k in 1:K){
		diag(SigmaINV[k,,]) <- 0.5 + 0.5*runif(p)
	}
	OmegaINV <- diag(q)
	# compute sufficient stats 
	suf_stat <- compute_sufficient_statistics_given_mu(y = y, 
	 z = z, K = K, x_data = x_data, mu = mu)

	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colMeans(x_data)
	priorConst1 <- array(diag(T_INV)*ksi, dim =c(p,1))
	# now simulate mu and Lambda
	f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_CUU(SigmaINV = SigmaINV, 
                suff_statistics = suf_stat, OmegaINV = OmegaINV, 
                K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
	# f2$mu contains the simulated means
	# f2$Lambdas contains the simulated factor loadings




cleanEx()
nameEx("compute_A_B_G_D_and_simulate_mu_Lambda_Sj")
### * compute_A_B_G_D_and_simulate_mu_Lambda_Sj

flush(stderr()); flush(stdout())

### Name: compute_A_B_G_D_and_simulate_mu_Lambda_Sj
### Title: Computation and simulations
### Aliases: compute_A_B_G_D_and_simulate_mu_Lambda_Sj

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- scale(as.matrix(waveDataset1500[ 1:20, -1])) # data
	z <-  waveDataset1500[ 1:20, 1] # class
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	q <- 2
	K <- length(table(z))           # 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	mu <- array( runif(K * p), dim = c(K,p) )
	y <- array(rnorm(n = q*n), dim = c(n,q))
	SigmaINV <- array(data = 0, dim = c(K,p,p) )
	for(k in 1:K){
		diag(SigmaINV[k,,]) <- 0.5 + 0.5*runif(p)
	}
	OmegaINV <- diag(q)
	# compute sufficient stats 
	suf_stat <- compute_sufficient_statistics(y = y, 
	 z = z, K = K, x_data = x_data)

	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colMeans(x_data)
	priorConst1 <- array(diag(T_INV)*ksi, dim =c(p,1))
	# now simulate mu and Lambda
	f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_Sj(SigmaINV = SigmaINV, 
                suff_statistics = suf_stat, OmegaINV = OmegaINV, 
                K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
	# f2$mu contains the simulated means
	# f2$Lambdas contains the simulated factor loadings




cleanEx()
nameEx("compute_sufficient_statistics")
### * compute_sufficient_statistics

flush(stderr()); flush(stdout())

### Name: compute_sufficient_statistics
### Title: Compute sufficient statistics
### Aliases: compute_sufficient_statistics

### ** Examples

        data(waveDataset1500)
        x_data <- as.matrix(waveDataset1500[ 1:20, -1]) # data
        z <-  waveDataset1500[ 1:20, 1] # class
        p <- dim(x_data)[2]
        n <- dim(x_data)[1]
        q <- 2
        K <- length(table(z))           # 3 classes
        # give some arbitrary values to the parameters:
        set.seed(1)
	y <- array(rnorm(n = q*n), dim = c(n,q))
	# compute sufficient stats 
	suf_stat <- compute_sufficient_statistics(y = y, 
	 z = z, K = K, x_data = x_data)



cleanEx()
nameEx("compute_sufficient_statistics_given_mu")
### * compute_sufficient_statistics_given_mu

flush(stderr()); flush(stdout())

### Name: compute_sufficient_statistics_given_mu
### Title: Compute sufficient statistics given mu
### Aliases: compute_sufficient_statistics_given_mu

### ** Examples

        data(waveDataset1500)
        x_data <- as.matrix(waveDataset1500[ 1:20, -1]) # data
        z <-  waveDataset1500[ 1:20, 1] # class
        p <- dim(x_data)[2]
        n <- dim(x_data)[1]
        q <- 2
        K <- length(table(z))           # 3 classes
        # give some arbitrary values to the parameters:
        set.seed(1)
        mu <- array( runif(K * p), dim = c(K,p) )
	y <- array(rnorm(n = q*n), dim = c(n,q))
	# compute sufficient stats 
	suf_stat <- compute_sufficient_statistics_given_mu(y = y, 
	 z = z, K = K, x_data = x_data, mu = mu)



cleanEx()
nameEx("compute_sufficient_statistics_q0")
### * compute_sufficient_statistics_q0

flush(stderr()); flush(stdout())

### Name: compute_sufficient_statistics_q0
### Title: Compute sufficient statistics for q = 0
### Aliases: compute_sufficient_statistics_q0

### ** Examples

        data(waveDataset1500)
        x_data <- as.matrix(waveDataset1500[ 1:20, -1]) # data
        z <-  waveDataset1500[ 1:20, 1] # class
        p <- dim(x_data)[2]
        n <- dim(x_data)[1]
        q <- 2
        K <- length(table(z))           # 3 classes
	# compute sufficient stats 
	suf_stat <- compute_sufficient_statistics_q0(
	 z = z, K = K, x_data = x_data)



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

# TOY EXAMPLE (very small numbers... only for CRAN check purposes)

#################################################################
# (a) using 2 cores in parallel, each one running 2 heated chains.
#################################################################
library('fabMix')

n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2		     # true number of clusters

sINV_diag = 1/((1:p))	 # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
			sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)

# Run `fabMix` for a _small_ number of iterations for the 
#	`UUU` (maximal model) and `CCC` (minimal model) parameterizations,
# 	using the default prior parallel heating parameters `dirPriorAlphas`.
#	NOTE: `dirPriorAlphas` may require some tuning in general.


qRange <- 2	# values for the number of factors (only the true number 
#                                                    is considered here)
Kmax <- 4	# number of components for the overfitted mixture model
nChains <- 2	# number of parallel heated chains

set.seed(1)
fm <- fabMix( model = c("UUU", "CCC"), nChains = nChains, 
	rawData = syntheticDataset$data, outDir = "toyExample",
        Kmax = Kmax, mCycles = 4, burnCycles = 1, q = qRange,
        g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
        warm_up_overfitting = 2, warm_up = 5) 

# WARNING: the following parameters: 
#  Kmax, nChains, mCycles, burnCycles, warm_up_overfitting, warm_up 
#	 should take (much) _larger_ values. E.g. a typical implementation consists of:
#        Kmax = 20, nChains >= 3, mCycles = 1100, burnCycles = 100, 
#        warm_up_overfitting = 500, warm_up = 5000. 

# Now print a run summary and produce some plots. 
print(fm)
plot(fm, what = "BIC")

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
##D # the next command takes ~ 1 hour in a Linux workstation with 12 threads.
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

# TOY EXAMPLE (very small numbers... only for CRAN check purposes)

#################################################################
# (a) using 2 cores in parallel, each one running 2 heated chains.
#################################################################
library('fabMix')

n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2		     # true number of clusters

sINV_diag = 1/((1:p))	 # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
			sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)

# Run `fabMix` for a _small_ number of iterations for the 
#	`UUU` (maximal model) and `CCC` (minimal model) parameterizations,
# 	using the default prior parallel heating parameters `dirPriorAlphas`.
#	NOTE: `dirPriorAlphas` may require some tuning in general.


qRange <- 2	# values for the number of factors (only the true number 
#                                                    is considered here)
Kmax <- 4	# number of components for the overfitted mixture model
nChains <- 2	# number of parallel heated chains

set.seed(1)
fm <- fabMix( model = c("UUU", "CCC"), nChains = nChains, 
	rawData = syntheticDataset$data, outDir = "toyExample",
        Kmax = Kmax, mCycles = 4, burnCycles = 1, q = qRange,
        g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
        warm_up_overfitting = 2, warm_up = 5) 

# WARNING: the following parameters: 
#  Kmax, nChains, mCycles, burnCycles, warm_up_overfitting, warm_up 
#	 should take (much) _larger_ values. E.g. a typical implementation consists of:
#        Kmax = 20, nChains >= 3, mCycles = 1100, burnCycles = 100, 
#        warm_up_overfitting = 500, warm_up = 5000. 

# You may also print and plot
# print(fm)
# plot(fm, what = "BIC")

#################################################################
# (b) using 12 cores_____________________________________________
#_______4 models with 3 heated chains running in parallel________
#_______considering all 8 model parameterizations________________
#################################################################
## Not run: 
##D library('fabMix')
##D set.seed(99)
##D n = 200                # sample size
##D p = 30                # number of variables
##D q = 2                # number of factors
##D K = 5		     # number of clusters
##D sINV_diag = rep(1/20,p) 	# diagonal of inverse variance of errors
##D syntheticDataset <- simData(sameLambda=FALSE,K.true = K, n = n, q = q, p = p, 
##D 			sINV_values = sINV_diag)
##D colnames(syntheticDataset$data) <- paste0("x_",1:p)
##D qRange <- 1:3	# range of values for the number of factors
##D Kmax <- 20	# number of components for the overfitted mixture model
##D nChains <- 3	# number of parallel heated chains
##D 
##D # the next command takes ~ 2 hours in a Linux machine with 12 threads.
##D 
##D fm <- fabMix( parallelModels = 4, 
##D 	nChains = nChains, 
##D 	model = c("UUU","CUU","UCU","CCU","UCC","UUC","CUC","CCC"), 
##D 	rawData = syntheticDataset$data, outDir = "toyExample_b",
##D         Kmax = Kmax, mCycles = 1100, burnCycles = 100, q = qRange) 
##D 
##D print(fm)
##D plot(fm, what = "BIC")
##D plot(fm, what = "classification_pairs")
##D # see also
##D # plot(fm); summary(fm)
##D 
## End(Not run)





cleanEx()
nameEx("observed.log.likelihood0")
### * observed.log.likelihood0

flush(stderr()); flush(stdout())

### Name: observed.log.likelihood0
### Title: Log-likelihood of the mixture model
### Aliases: observed.log.likelihood0

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- waveDataset1500[ 1:20, -1] # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( runif(K * p), dim = c(K,p) )
	Lambda <- array( runif(K*p*q), dim = c(K,p,q) )
	SigmaINV <- array(1, dim = c(p,p))
	Sigma <- 1/diag(SigmaINV)
	# compute the complete.log.likelihood
	observed.log.likelihood0(x_data = x_data, w = w, 
		mu = mu, Lambda = Lambda, Sigma = Sigma, z = z)



cleanEx()
nameEx("observed.log.likelihood0_Sj")
### * observed.log.likelihood0_Sj

flush(stderr()); flush(stdout())

### Name: observed.log.likelihood0_Sj
### Title: Log-likelihood of the mixture model
### Aliases: observed.log.likelihood0_Sj

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- waveDataset1500[ 1:20, -1] # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( runif(K * p), dim = c(K,p) )
	Lambda <- array( runif(K*p*q), dim = c(K,p,q) )
	Sigma <- matrix(1:K, nrow = K, ncol = p)
	# compute the complete.log.likelihood
	observed.log.likelihood0_Sj(x_data = x_data, w = w, 
		mu = mu, Lambda = Lambda, Sigma = Sigma, z = z)



cleanEx()
nameEx("observed.log.likelihood0_Sj_q0")
### * observed.log.likelihood0_Sj_q0

flush(stderr()); flush(stdout())

### Name: observed.log.likelihood0_Sj_q0
### Title: Log-likelihood of the mixture model for q=0
### Aliases: observed.log.likelihood0_Sj_q0

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- waveDataset1500[ 1:20, -1] # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( runif(K * p), dim = c(K,p) )
	Sigma <- matrix(1:K, nrow = K, ncol = p)
	# compute the complete.log.likelihood
	observed.log.likelihood0_Sj_q0(x_data = x_data, w = w, 
		mu = mu, Sigma = Sigma, z = z)



cleanEx()
nameEx("observed.log.likelihood0_q0_sameSigma")
### * observed.log.likelihood0_q0_sameSigma

flush(stderr()); flush(stdout())

### Name: observed.log.likelihood0_q0_sameSigma
### Title: Log-likelihood of the mixture model for q=0 and same variance of
###   errors
### Aliases: observed.log.likelihood0_q0_sameSigma

### ** Examples

	library('fabMix')
	data(waveDataset1500)
	x_data <- waveDataset1500[ 1:20, -1] # data
	z <-  waveDataset1500[ 1:20, 1]	# class
	p <- dim(x_data)[2]
	q <- 2
	K <- length(table(z))		# 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	w <- rep(1/K, K)
	mu <- array( runif(K * p), dim = c(K,p) )
	SigmaINV <- array(1, dim = c(p,p))
	Sigma <- 1/diag(SigmaINV)
	# compute the complete.log.likelihood
	observed.log.likelihood0_q0_sameSigma(x_data = x_data, w = w, 
		mu = mu,  Sigma = Sigma, z = z)



cleanEx()
nameEx("overfittingMFA")
### * overfittingMFA

flush(stderr()); flush(stdout())

### Name: overfittingMFA
### Title: Basic MCMC sampler for the 'UCU' model
### Aliases: overfittingMFA

### ** Examples

library('fabMix')
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
Kmax <- 4       # number of components for the overfitted mixture model

set.seed(1)
overfittingMFA(x_data = syntheticDataset$data, 
	originalX = syntheticDataset$data, outputDirectory = 'outDir', 
	Kmax = Kmax, m = 5, burn = 1, 
	g = 0.5, h = 0.5, alpha_prior = rep(1, Kmax), 
	alpha_sigma = 0.5, beta_sigma = 0.5, 
	start_values = FALSE, q = 2,  gibbs_z = 1)
list.files('outDir')
unlink('outDir', recursive = TRUE)




cleanEx()
nameEx("overfittingMFA_CCC")
### * overfittingMFA_CCC

flush(stderr()); flush(stdout())

### Name: overfittingMFA_CCC
### Title: Basic MCMC sampler for the 'CCC' model
### Aliases: overfittingMFA_CCC

### ** Examples

library('fabMix')
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
Kmax <- 4       # number of components for the overfitted mixture model

set.seed(1)
overfittingMFA_CCC <- overfittingMFA_CCC(x_data = syntheticDataset$data, 
	originalX = syntheticDataset$data, outputDirectory = 'outDir', 
	Kmax = Kmax, m = 5, burn = 1, 
	g = 0.5, h = 0.5, alpha_prior = rep(1, Kmax), 
	alpha_sigma = 0.5, beta_sigma = 0.5, 
	start_values = FALSE, q = 2,  gibbs_z = 1)
list.files('outDir')
unlink('outDir', recursive = TRUE)




cleanEx()
nameEx("overfittingMFA_CCU")
### * overfittingMFA_CCU

flush(stderr()); flush(stdout())

### Name: overfittingMFA_CCU
### Title: Basic MCMC sampler for the 'CCU' model
### Aliases: overfittingMFA_CCU

### ** Examples

library('fabMix')
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
Kmax <- 4       # number of components for the overfitted mixture model

set.seed(1)
overfittingMFA_CCU(x_data = syntheticDataset$data, 
	originalX = syntheticDataset$data, outputDirectory = 'outDir', 
	Kmax = Kmax, m = 5, burn = 1, 
	g = 0.5, h = 0.5, alpha_prior = rep(1, Kmax), 
	alpha_sigma = 0.5, beta_sigma = 0.5, 
	start_values = FALSE, q = 2,  gibbs_z = 1)
list.files('outDir')
unlink('outDir', recursive = TRUE)




cleanEx()
nameEx("overfittingMFA_CUC")
### * overfittingMFA_CUC

flush(stderr()); flush(stdout())

### Name: overfittingMFA_CUC
### Title: Basic MCMC sampler for the 'CUC' model
### Aliases: overfittingMFA_CUC

### ** Examples

library('fabMix')
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
Kmax <- 4       # number of components for the overfitted mixture model

set.seed(1)
overfittingMFA_CUC(x_data = syntheticDataset$data, 
	originalX = syntheticDataset$data, outputDirectory = 'outDir', 
	Kmax = Kmax, m = 5, burn = 1, 
	g = 0.5, h = 0.5, alpha_prior = rep(1, Kmax), 
	alpha_sigma = 0.5, beta_sigma = 0.5, 
	start_values = FALSE, q = 2,  gibbs_z = 1)
list.files('outDir')
unlink('outDir', recursive = TRUE)




cleanEx()
nameEx("overfittingMFA_CUU")
### * overfittingMFA_CUU

flush(stderr()); flush(stdout())

### Name: overfittingMFA_CUU
### Title: Basic MCMC sampler for the 'CUU' model
### Aliases: overfittingMFA_CUU

### ** Examples

library('fabMix')
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
Kmax <- 4       # number of components for the overfitted mixture model

set.seed(1)
overfittingMFA_CUU(x_data = syntheticDataset$data, 
	originalX = syntheticDataset$data, outputDirectory = 'outDir', 
	Kmax = Kmax, m = 5, burn = 1, 
	g = 0.5, h = 0.5, alpha_prior = rep(1, Kmax), 
	alpha_sigma = 0.5, beta_sigma = 0.5, 
	start_values = FALSE, q = 2,  gibbs_z = 1)
list.files('outDir')
unlink('outDir', recursive = TRUE)




cleanEx()
nameEx("overfittingMFA_Sj")
### * overfittingMFA_Sj

flush(stderr()); flush(stdout())

### Name: overfittingMFA_Sj
### Title: Basic MCMC sampler for the 'UUU' model
### Aliases: overfittingMFA_Sj

### ** Examples

library('fabMix')
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
Kmax <- 4       # number of components for the overfitted mixture model

set.seed(1)
overfittingMFA_Sj(x_data = syntheticDataset$data, 
	originalX = syntheticDataset$data, outputDirectory = 'outDir', 
	Kmax = Kmax, m = 5, burn = 1, 
	g = 0.5, h = 0.5, alpha_prior = rep(1, Kmax), 
	alpha_sigma = 0.5, beta_sigma = 0.5, 
	start_values = FALSE, q = 2,  gibbs_z = 1)
list.files('outDir')
unlink('outDir', recursive = TRUE)




cleanEx()
nameEx("overfittingMFA_UCC")
### * overfittingMFA_UCC

flush(stderr()); flush(stdout())

### Name: overfittingMFA_UCC
### Title: Basic MCMC sampler for the 'UCC' model
### Aliases: overfittingMFA_UCC

### ** Examples

library('fabMix')
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
Kmax <- 4       # number of components for the overfitted mixture model

set.seed(1)
overfittingMFA_UCC(x_data = syntheticDataset$data, 
	originalX = syntheticDataset$data, outputDirectory = 'outDir', 
	Kmax = Kmax, m = 5, burn = 1, 
	g = 0.5, h = 0.5, alpha_prior = rep(1, Kmax), 
	alpha_sigma = 0.5, beta_sigma = 0.5, 
	start_values = FALSE, q = 2,  gibbs_z = 1)
list.files('outDir')
unlink('outDir', recursive = TRUE)




cleanEx()
nameEx("overfittingMFA_UUC")
### * overfittingMFA_UUC

flush(stderr()); flush(stdout())

### Name: overfittingMFA_UUC
### Title: Basic MCMC sampler for the 'UUC' model
### Aliases: overfittingMFA_UUC

### ** Examples

library('fabMix')
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
colnames(syntheticDataset$data) <- paste0("x_",1:p)
Kmax <- 4       # number of components for the overfitted mixture model

set.seed(1)
overfittingMFA_UUC(x_data = syntheticDataset$data, 
	originalX = syntheticDataset$data, outputDirectory = 'outDir', 
	Kmax = Kmax, m = 5, burn = 1, 
	g = 0.5, h = 0.5, alpha_prior = rep(1, Kmax), 
	alpha_sigma = 0.5, beta_sigma = 0.5, 
	start_values = FALSE, q = 2,  gibbs_z = 1)
list.files('outDir')
unlink('outDir', recursive = TRUE)




cleanEx()
nameEx("simData")
### * simData

flush(stderr()); flush(stdout())

### Name: simData
### Title: Synthetic data generator
### Aliases: simData

### ** Examples

library('fabMix')

n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
summary(syntheticDataset)



cleanEx()
nameEx("simData2")
### * simData2

flush(stderr()); flush(stdout())

### Name: simData2
### Title: Synthetic data generator 2
### Aliases: simData2

### ** Examples

library('fabMix')

n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData2(K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
summary(syntheticDataset)



cleanEx()
nameEx("update_OmegaINV")
### * update_OmegaINV

flush(stderr()); flush(stdout())

### Name: update_OmegaINV
### Title: Gibbs sampling for Omega^{-1}
### Aliases: update_OmegaINV

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
SigmaINV <- array(data = 0, dim = c(K,p,p))
for(k in 1:K){
        diag(SigmaINV[k,,]) <- 1/diag(syntheticDataset$variance) + rgamma(p, shape=1, rate = 1)
}

# use the real values as input and simulate allocations
update_OmegaINV(Lambda = syntheticDataset$factorLoadings, 
        K = K, g=0.5, h = 0.5)




cleanEx()
nameEx("update_OmegaINV_Cxx")
### * update_OmegaINV_Cxx

flush(stderr()); flush(stdout())

### Name: update_OmegaINV_Cxx
### Title: Gibbs sampling for Omega^{-1} for Cxx model
### Aliases: update_OmegaINV_Cxx

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
SigmaINV <- array(data = 0, dim = c(K,p,p))
for(k in 1:K){
        diag(SigmaINV[k,,]) <- 1/diag(syntheticDataset$variance) + rgamma(p, shape=1, rate = 1)
}

# Use the real values as input and simulate allocations.
# Mmake sure that in this case Lambda[k,,] is the same  
# for all k = 1,..., K
update_OmegaINV_Cxx(Lambda = syntheticDataset$factorLoadings, 
        K = K, g=0.5, h = 0.5)




cleanEx()
nameEx("update_SigmaINV_faster")
### * update_SigmaINV_faster

flush(stderr()); flush(stdout())

### Name: update_SigmaINV_faster
### Title: Gibbs sampling for Sigma^{-1}
### Aliases: update_SigmaINV_faster

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)

# use the real values as input and update SigmaINV
update_SigmaINV_faster(x_data = syntheticDataset$data, 
	z = syntheticDataset$class, 
	y = syntheticDataset$factors, 
	Lambda = syntheticDataset$factorLoadings, 
	mu = syntheticDataset$means, 
	K = K, 
	alpha_sigma = 0.5, beta_sigma = 0.5)




cleanEx()
nameEx("update_SigmaINV_faster_Sj")
### * update_SigmaINV_faster_Sj

flush(stderr()); flush(stdout())

### Name: update_SigmaINV_faster_Sj
### Title: Gibbs sampling for Sigma^{-1} per component
### Aliases: update_SigmaINV_faster_Sj

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)

# use the real values as input and update SigmaINV
update_SigmaINV_faster_Sj(x_data = syntheticDataset$data, 
	z = syntheticDataset$class, 
	y = syntheticDataset$factors, 
	Lambda = syntheticDataset$factorLoadings, 
	mu = syntheticDataset$means, 
	K = K, 
	alpha_sigma = 0.5, beta_sigma = 0.5)




cleanEx()
nameEx("update_SigmaINV_xCC")
### * update_SigmaINV_xCC

flush(stderr()); flush(stdout())

### Name: update_SigmaINV_xCC
### Title: Gibbs sampling for Sigma^{-1} for xCC models
### Aliases: update_SigmaINV_xCC

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)

# use the real values as input and update SigmaINV
update_SigmaINV_xCC(x_data = syntheticDataset$data, 
	z = syntheticDataset$class, 
	y = syntheticDataset$factors, 
	Lambda = syntheticDataset$factorLoadings, 
	mu = syntheticDataset$means, 
	K = K, 
	alpha_sigma = 0.5, beta_sigma = 0.5)




cleanEx()
nameEx("update_SigmaINV_xUC")
### * update_SigmaINV_xUC

flush(stderr()); flush(stdout())

### Name: update_SigmaINV_xUC
### Title: Gibbs sampling for Sigma^{-1} per component for xUC models
### Aliases: update_SigmaINV_xUC

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)

# use the real values as input and update SigmaINV
update_SigmaINV_xUC(x_data = syntheticDataset$data, 
	z = syntheticDataset$class, 
	y = syntheticDataset$factors, 
	Lambda = syntheticDataset$factorLoadings, 
	mu = syntheticDataset$means, 
	K = K, 
	alpha_sigma = 0.5, beta_sigma = 0.5)




cleanEx()
nameEx("update_all_y")
### * update_all_y

flush(stderr()); flush(stdout())

### Name: update_all_y
### Title: Gibbs sampling for y in 'xCx' model
### Aliases: update_all_y

### ** Examples

library('fabMix')

n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
# use the real values as input and simulate factors
update_all_y(x_data = syntheticDataset$data, 
		mu = syntheticDataset$means, 
		SigmaINV = diag(1/diag(syntheticDataset$variance)), 
		Lambda = syntheticDataset$factorLoadings, 
		z = syntheticDataset$class)




cleanEx()
nameEx("update_all_y_Sj")
### * update_all_y_Sj

flush(stderr()); flush(stdout())

### Name: update_all_y_Sj
### Title: Gibbs sampling for y in 'xUx' model
### Aliases: update_all_y_Sj

### ** Examples

library('fabMix')

n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
# add some noise here:
SigmaINV <- array(data = 0, dim = c(K,p,p))
for(k in 1:K){
        diag(SigmaINV[k,,]) <- 1/diag(syntheticDataset$variance) + rgamma(p, shape=1, rate = 1)
}

# use the real values as input and simulate factors
update_all_y_Sj(x_data = syntheticDataset$data, 
		mu = syntheticDataset$means, 
		SigmaINV = SigmaINV, 
		Lambda = syntheticDataset$factorLoadings, 
		z = syntheticDataset$class)




cleanEx()
nameEx("update_z2")
### * update_z2

flush(stderr()); flush(stdout())

### Name: update_z2
### Title: Collapsed Gibbs for z using matrix inversion lemma
### Aliases: update_z2

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
# use the real values as input and simulate allocations
update_z2(w = syntheticDataset$weights, mu = syntheticDataset$means, 
	Lambda = syntheticDataset$factorLoadings, 
	SigmaINV = diag(1/diag(syntheticDataset$variance)), 
	K = K, x_data = syntheticDataset$data)$z




cleanEx()
nameEx("update_z2_Sj")
### * update_z2_Sj

flush(stderr()); flush(stdout())

### Name: update_z2_Sj
### Title: Collapsed Gibbs for z using matrix inversion lemma
### Aliases: update_z2_Sj

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
SigmaINV <- array(data = 0, dim = c(K,p,p))
for(k in 1:K){
	diag(SigmaINV[k,,]) <- 1/diag(syntheticDataset$variance) + rgamma(p, shape=1, rate = 1)
}

# use the real values as input and simulate allocations
update_z2_Sj(w = syntheticDataset$weights, mu = syntheticDataset$means, 
	Lambda = syntheticDataset$factorLoadings, 
	SigmaINV = SigmaINV, 
	K = K, x_data = syntheticDataset$data)$z




cleanEx()
nameEx("update_z4")
### * update_z4

flush(stderr()); flush(stdout())

### Name: update_z4
### Title: Collapsed Gibbs for z
### Aliases: update_z4

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
# use the real values as input and simulate allocations
update_z4(w = syntheticDataset$weights, mu = syntheticDataset$means, 
	Lambda = syntheticDataset$factorLoadings, 
	SigmaINV = diag(1/diag(syntheticDataset$variance)), 
	K = K, x_data = syntheticDataset$data)$z




cleanEx()
nameEx("update_z4_Sj")
### * update_z4_Sj

flush(stderr()); flush(stdout())

### Name: update_z4_Sj
### Title: Collapsed Gibbs for z
### Aliases: update_z4_Sj

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
SigmaINV <- array(data = 0, dim = c(K,p,p))
for(k in 1:K){
	diag(SigmaINV[k,,]) <- 1/diag(syntheticDataset$variance) + rgamma(p, shape=1, rate = 1)
}

# use the real values as input and simulate allocations
update_z4_Sj(w = syntheticDataset$weights, mu = syntheticDataset$means, 
	Lambda = syntheticDataset$factorLoadings, 
	SigmaINV = SigmaINV, 
	K = K, x_data = syntheticDataset$data)$z




cleanEx()
nameEx("update_z_b")
### * update_z_b

flush(stderr()); flush(stdout())

### Name: update_z_b
### Title: Gibbs sampling for z
### Aliases: update_z_b

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)

# use the real values as input and simulate allocations
update_z_b(w = syntheticDataset$weights, mu = syntheticDataset$means, 
	Lambda = syntheticDataset$factorLoadings, 
	y = syntheticDataset$factors,
	SigmaINV = diag(1/diag(syntheticDataset$variance)), 
	K = K, x_data = syntheticDataset$data)$z




cleanEx()
nameEx("update_z_b_Sj")
### * update_z_b_Sj

flush(stderr()); flush(stdout())

### Name: update_z_b_Sj
### Title: Gibbs sampling for z
### Aliases: update_z_b_Sj

### ** Examples

library('fabMix')
# simulate some data
n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters
sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
SigmaINV <- array(data = 0, dim = c(K,p,p))
for(k in 1:K){
	diag(SigmaINV[k,,]) <- 1/diag(syntheticDataset$variance) + rgamma(p, shape=1, rate = 1)
}

# use the real values as input and simulate allocations
update_z_b_Sj(w = syntheticDataset$weights, mu = syntheticDataset$means, 
	Lambda = syntheticDataset$factorLoadings, 
	y = syntheticDataset$factors,
	SigmaINV = SigmaINV, 
	K = K, x_data = syntheticDataset$data)$z




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
