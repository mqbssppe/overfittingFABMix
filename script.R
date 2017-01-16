
Kmax <- 20
#alphas <- c(40, 30, 20, 10, 5, 3, 1, 0.5, 0.5^2, 0.5^3, 0.5^4, 0.05)
set.seed(10)
#source('simMix.R')
source('sim2.R')
#c(1,10, 20, 30, 100)/Kmax




alphas <- c(c(20, 10, 5, 1, 0.50, 0.25),c((1/Kmax)*seq(3,1.05, length = 6)), c((1/Kmax)*seq(1.04,1, length = 4)))
alphas <- alphas[length(alphas):1]

source('~/Dropbox/sparseFA_MIX/heated_prior/csf_tdp.R')
#	q <- 2	
for( q in 1:10 ){
	p <- dim(originalX)[2]
	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}
	d = 2*p + p*q + q*(q-1)/2
	heated_chains(dirPriorAlphas = alphas, rawData = originalX, outDir = paste0('q_',q), 
			Kmax = Kmax, mCycles = 1050, burnCycles = 10, g = 0.2, h = 0.2, alpha_sigma = 5, beta_sigma = 1, q = q) 

}
## ftiakse ta arguments 

source('~/Dropbox/sparseFA_MIX/conjugate_not_sparse.R')

truncatedDirichletProcess(outputDirectory = 'test_conjugeate', Kmax = Kmax, m = 1000, thinning = 10, burn = 10, c_1 = 0.1, c_2 = 0.5, alpha_prior = rep(1/Kmax, Kmax), alpha_sigma = 0.5, beta_sigma = 0.5, progressGraphs = T)




library('mixAK')

mixakModels2 <- vector('list', length = 20)
for (k in 1:20){
 Prior1 <- list(priorK = "fixed", Kmax = k, delta = 1, priormuQ = "independentC")
 mixakModels2[[k]] <- NMixMCMC(y0 = x_data, prior = Prior1, nMCMC = c(burn = 5000, keep = 25000, thin = 10, info = 1000),  PED = TRUE, parallel = TRUE)
}



alphas <- c(seq(1, 0.25, length = 4),c((1/Kmax)*seq(3,1.05, length = 4)), c((1/Kmax)*seq(1.04,1, length = 2)))
alphas <- alphas[length(alphas):1]

source('~/Dropbox/sparseFA_MIX/heated_prior/csf_tdp_same_sigma.R')
#	q <- 2	
for( q in 1:10 ){
	p <- dim(originalX)[2]
	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}
	Kmax <- 20
#	dirPriorAlphas <- c(seq(1, 0.25, length = 4),c((1/Kmax)*seq(3,1.05, length = 4)), c((1/Kmax)*seq(1.04,1, length = 2)))
#	dirPriorAlphas <- dirPriorAlphas[length(dirPriorAlphas):1]
#	dirPriorAlphas <- dirPriorAlphas/3

	heated_chains( Kmax = 50, rawData = originalX, outDir = paste0('q_',q), mCycles = 5050, burnCycles = 10, g = 0.2, h = 0.2, alpha_sigma = 1, beta_sigma = 1, q = q) 

}





	meanR <- mean(diag(var(x_data)))
	p <- dim(originalX)[2]
	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}
	heated_chains(dirPriorAlphas = alphas, rawData = originalX, outDir = paste0('q_',q), 
			Kmax = Kmax, mCycles = 25, burnCycles = 2, g = 0.2, h = 0.2, alpha_sigma = 1, beta_sigma = 1, q = q) 

	meanR <- mean(diag(var(originalX)))
	heated_chains(dirPriorAlphas = alphas, rawData = originalX, outDir = paste0('q_',q), 
			Kmax = Kmax, mCycles = 550, burnCycles = 10, g = 0.2, h = 0.2*meanR, alpha_sigma = 1, beta_sigma =  1*meanR, q = q, normalize = FALSE) 


source('~/Dropbox/sparseFA_MIX/heated_prior/q_0.R')
heated_chains(dirPriorAlphas = alphas, rawData = originalX, outDir = paste0('q_',q),Kmax = Kmax, mCycles = 1000, burnCycles = 20, g = 1, h = 1, alpha_sigma = 1, beta_sigma = 1, q = q)
Km <- 20
outputDirectory <- paste0('q_',q)
source('~/Dropbox/sparseFA_MIX/heated_prior/dic_full_q_0.R')

source('~/Dropbox/sparseFA_MIX/heated_prior/q_0_same_sigma.R')
heated_chains(dirPriorAlphas = alphas, rawData = originalX, outDir = paste0('same_q_',q),Kmax = Kmax, mCycles = 100, burnCycles = 20, g = 0.2, h = 0.2, alpha_sigma = 1, beta_sigma = 1, q = q)
Km <- 20
outputDirectory <- paste0('same_q_',q)
source('~/Dropbox/sparseFA_MIX/heated_prior/dic_full_same_sigma_q_0.R')

z <- as.matrix(read.table("q_4/zValues.txt")); ari <- apply(z,1,function(y){adjustedRandIndex(z.true,y)});plot(ari)
abline(h = adjustedRandIndex(z.true, mod2$classification))


	source('~/Dropbox/sparseFA_MIX/heated_prior/csf_tdp_same_sigma.R')
	q <- 4
	p <- dim(originalX)[2]
	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}
	Kmax <- 20
	alphas <- (1/Kmax)*exp(seq(0,2,length = 10))
	alphas <- alphas[1:8]
#	dirPriorAlphas <- dirPriorAlphas[length(dirPriorAlphas):1]
#	dirPriorAlphas <- dirPriorAlphas/3
#	alphas <- c(1/Kmax, 1)	
	heated_chains( dirPriorAlphas = alphas, rawData = originalX, outDir = paste0('q_',q), mCycles = 5050, burnCycles = 10, g = 0.5, h = 0.5, alpha_sigma = 0.1, beta_sigma = 0.1, q = q) 
heated_chains( rawData = originalX, outDir = paste0('q_',q), mCycles = 5050, burnCycles = 10, g = 1, h = 1, alpha_sigma = 0.2, beta_sigma = 0.2, q = q) 
#	heated_chains( rawData = originalX, outDir = paste0('q_',q), mCycles = 5050, burnCycles = 10, g = 0.2, h = 0.2, alpha_sigma = 0.2, beta_sigma = 0.2, q = q) 
