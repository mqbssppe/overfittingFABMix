library(fabMix)
library(pgmm)
library(EMMIXskew)
library(EMMIXmfa)
library(flexmix)
library(mclust)


n = 1000              # sample size
p = 20                # number of variables
q = 4                 # number of factors
K = 2                # number of clusters
#			scenario 1:
sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
#			or you can also use scenario 2:
#sINV_diag = rgamma(p, shape = 1, rate = 1)
#			scenario 1 is more challenging
#			as it produces more noisy observations
#			for features j --> p 
set.seed(10)
syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )


Kmax <- 20

#	rG: range of mixture components
#	rq: range of number of factors
#	zStart=2 denotes kMeans initialization
#	loop is only relevant when zStart = 1	
pgmmRun <- pgmmEM(x = syntheticDataset$data, rG=1:Kmax, rq=1:5, class=NULL, icl=FALSE, zstart=2, 
	cccStart=TRUE, loop=1, zlist=NULL,
	modelSubset=c("UCU"), tol=0.1, relax=FALSE)

pgmmK <- pgmmRun$g
pgmmZ <- apply(pgmmRun$zhat, 1, function(y){order(y, decreasing=TRUE)[1]})
adjustedRandIndex(syntheticDataset$class, apply(pgmmRun$zhat, 1, function(y){order(y, decreasing=TRUE)[1]}))

# Em Mix einai cool (default covariance corresponds to ncov = 3 with general cov structure)

EmSkewRun1 <- EmSkewRun2 <- EmSkewRun3 <- EmSkewRun4 <- EmSkewRun5 <- vector("list", length = Kmax)
for(k in 1:Kmax){
	cat(paste0("EmMix run with k = ", k),"\n")
	cat(paste0("      ncov = ", 1),"\n")
	EmSkewRun1[[k]] <- EmSkew(dat = syntheticDataset$data, g = k, distr="mvn", ncov=1,clust=NULL,itmax=2000,
			epsilon=1e-6, nkmeans=0, nrandom=10,nhclust=FALSE,debug=FALSE,initloop=20)
	cat(paste0("      ncov = ", 2),"\n")
	EmSkewRun2[[k]] <- EmSkew(dat = syntheticDataset$data, g = k, distr="mvn", ncov=2,clust=NULL,itmax=2000,
			epsilon=1e-6, nkmeans=0, nrandom=10,nhclust=FALSE,debug=FALSE,initloop=20)
	cat(paste0("      ncov = ", 3),"\n")
	EmSkewRun3[[k]] <- EmSkew(dat = syntheticDataset$data, g = k, distr="mvn", ncov=3,clust=NULL,itmax=2000,
			epsilon=1e-6, nkmeans=0, nrandom=10,nhclust=FALSE,debug=FALSE,initloop=20)
	cat(paste0("      ncov = ", 4),"\n")
	EmSkewRun4[[k]] <- EmSkew(dat = syntheticDataset$data, g = k, distr="mvn", ncov=4,clust=NULL,itmax=2000,
			epsilon=1e-6, nkmeans=0, nrandom=10,nhclust=FALSE,debug=FALSE,initloop=20)
	cat(paste0("      ncov = ", 5),"\n")
	EmSkewRun5[[k]] <- EmSkew(dat = syntheticDataset$data, g = k, distr="mvn", ncov=5,clust=NULL,itmax=2000,
			epsilon=1e-6, nkmeans=0, nrandom=10,nhclust=FALSE,debug=FALSE,initloop=20)
}
bicValues <- array(data = NA , dim = c(Kmax, 5))
bicValues[,1] <- unlist(lapply(EmSkewRun1, function(y){y$bic}))
bicValues[,2] <- unlist(lapply(EmSkewRun2, function(y){y$bic}))
bicValues[,3] <- unlist(lapply(EmSkewRun3, function(y){y$bic}))
bicValues[,4] <- unlist(lapply(EmSkewRun4, function(y){y$bic}))
bicValues[,5] <- unlist(lapply(EmSkewRun5, function(y){y$bic}))
matplot(bicValues[,-1], type = "l")
covModel <- apply(bicValues, 2, function(y){order(y, decreasing=F)[1]})
BestModels <- c(bicValues[covModel[1],1], bicValues[covModel[2],2], bicValues[covModel[3],3], bicValues[covModel[4],4], bicValues[covModel[5],5])
EmSkewCovModel <- order(BestModels, decreasing = FALSE)[1]
EmSkewK <- covModel[EmSkewCovModel]
AllModels <- vector("list", length = 5)
AllModels[[1]] <- EmSkewRun1
AllModels[[2]] <- EmSkewRun2
AllModels[[3]] <- EmSkewRun3
AllModels[[4]] <- EmSkewRun4
AllModels[[5]] <- EmSkewRun5
EmSkewZ <- AllModels[[EmSkewCovModel]][[EmSkewK]]$clust
adjustedRandIndex(EmSkewZ, syntheticDataset$class)


# run EMMIXmfa
EmMfaRun2 <- EmMfaRun3 <- EmMfaRun4 <- EmMfaRun5 <- EmMfaRun6 <- EmMfaRun7 <-  vector("list", length = Kmax)
for(k in 1:Kmax){
	EmMfaRun2[[k]] <- mfa(Y = syntheticDataset$data, g=k, q=2, itmax=1000, nkmeans=1, nrandom=5, sigmaType = "unique", Dtype = "common")
	EmMfaRun3[[k]] <- mfa(Y = syntheticDataset$data, g=k, q=3, itmax=1000, nkmeans=1, nrandom=5, sigmaType = "unique", Dtype = "common")
	EmMfaRun4[[k]] <- mfa(Y = syntheticDataset$data, g=k, q=4, itmax=1000, nkmeans=1, nrandom=5, sigmaType = "unique", Dtype = "common")
	EmMfaRun5[[k]] <- mfa(Y = syntheticDataset$data, g=k, q=5, itmax=1000, nkmeans=1, nrandom=5, sigmaType = "unique", Dtype = "common")
	EmMfaRun6[[k]] <- mfa(Y = syntheticDataset$data, g=k, q=6, itmax=1000, nkmeans=1, nrandom=5, sigmaType = "unique", Dtype = "common")
	EmMfaRun7[[k]] <- mfa(Y = syntheticDataset$data, g=k, q=7, itmax=1000, nkmeans=1, nrandom=5, sigmaType = "unique", Dtype = "common")	
}
bicValues <- array(data = NA , dim = c(Kmax, 6))
bicValues[,1] <- unlist(lapply(EmMfaRun2, function(y){y$BIC}))
bicValues[,2] <- unlist(lapply(EmMfaRun3, function(y){y$BIC}))
bicValues[,3] <- unlist(lapply(EmMfaRun4, function(y){y$BIC}))
bicValues[,4] <- unlist(lapply(EmMfaRun5, function(y){y$BIC}))
bicValues[,5] <- unlist(lapply(EmMfaRun6, function(y){y$BIC}))
bicValues[,6] <- unlist(lapply(EmMfaRun7, function(y){y$BIC}))
covModel <- apply(bicValues, 2, function(y){order(y, decreasing=F)[1]})
BestModels <- c(bicValues[covModel[1],1], bicValues[covModel[2],2], bicValues[covModel[3],3], bicValues[covModel[4],4], bicValues[covModel[5],5], bicValues[covModel[6],6])
EmMfaqModel <- order(BestModels, decreasing = FALSE)[1]
EmMfaK <- covModel[EmMfaqModel]
AllModels <- vector("list", length = 6)
AllModels[[1]] <- EmMfaRun2
AllModels[[2]] <- EmMfaRun3
AllModels[[3]] <- EmMfaRun4
AllModels[[4]] <- EmMfaRun5
AllModels[[5]] <- EmMfaRun6
AllModels[[6]] <- EmMfaRun7
EmMfaZ <- AllModels[[EmMfaqModel]][[EmMfaK]]$clust
#adjustedRandIndex(EmSkewZ, syntheticDataset$class)



