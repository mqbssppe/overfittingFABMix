library('label.switching')
#Km <- 20
dealWithLabelSwitching_S_j <- function(outputFolder, q, burn, z.true, Km){
	if(missing(Km)){Km <- 20}
	if(missing(burn)){burn <- 1:1}
	cat(paste0('-    (5) Dealing with label switching for q = ', q), '\n')
	setwd(outputFolder)
	cat(paste0('         * Entering directory: ', getwd()),'\n')
	z <- as.matrix(read.table("zValues.txt"))
	logl <- read.table("kValues.txt", header=T)
	logl <- logl[-burn, ]
	z <- z[-burn,]
#	print(table(logl[,1]))
	K <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
	kSelected <- K
	cat(paste0('         * Posterior mode corresponds to K = ', K),'\n')
	if(K == 1){
		cat(paste0('         *  no label switching algorithms are applied'),'\n')
	}else{
		index <- which(logl[,1] == K)
		Kindex <- index
		logl <- logl[index,2]
		z <- z[index,]
		mapIndex <- which(logl == max(logl))[1]
		zPivot <- as.numeric(z[mapIndex,])
		K <- max(z)
		if(K < 20){cat(paste('WARNING: max Z not equal to 20'),'\n')}
		if(K < Km){ 
			newZ <- t(apply(z, 1, function(y){ myI <- which(y == K); y[myI] <- rep(Km, length(myI));return(y) } ) )
			z <- newZ
			K <- max(z)
		}

		m <- length(logl)
		if(missing(z.true)){
			ls <- label.switching(method = c("ECR", "ECR-ITERATIVE-1"), zpivot = zPivot, z = z, K = K)}else{
			ls <- label.switching(method = c("ECR", "ECR-ITERATIVE-1"), zpivot = zPivot, z = z, K = K, groundTruth=z.true)
		}
		oldLabels <- 1:K
		index <- Kindex
		allocationsECR <- z
		for (i in 1:m){
			myPerm <- order(ls$permutations$"ECR"[i,])
			allocationsECR[i,] <- myPerm[z[i,]]
		}
		write.table(allocationsECR, file = 'reordered_allocations_ecr.txt')
		cat(paste0('         * write file: `reordered_allocations_ecr.txt`' ),'\n')
		l <- as.matrix(read.table(paste0("LambdaValues",1,".txt")))
		J <- dim(l)[2]
		mcmc <- array(data = NA, dim = c(m,K,J))
		for(k in 1:K){
		#       lMean <- apply(l,2,mean)
		#       lMean.matrix <- matrix(lMean,nrow = p, ncol = q, byrow=TRUE)
			l <- as.matrix(read.table(paste0("LambdaValues",k,".txt")))
			l <- l[-burn,]
			mcmc[,k,] <- l[index,]
		}
		lambda.perm.mcmc <- permute.mcmc(mcmc, ls$permutations$ECR)$output
		write.table(lambda.perm.mcmc, file = 'reordered_lambda_ecr.txt')
		cat(paste0('         * write file: `reordered_lambda_ecr.txt`'),'\n')
		lambda.mean <- array(data = NA, dim = c(K,p,q))
		lambda.map <- array(data = NA, dim = c(K,p,q))
		for(k in 1:K){
			lMean <- apply(lambda.perm.mcmc[,k,],2,mean)
			lambda.mean[k,,] <- matrix(lMean,nrow = p, ncol = q, byrow=TRUE)
			lambda.map[k,,] <- matrix(lambda.perm.mcmc[mapIndex,k,],nrow = p, ncol = q, byrow=TRUE)
		}
		write.table(lambda.mean, file = 'lambda_estimate_ecr.txt', col.names = paste('lambda',1:p, rep(1:q, each = p), sep = "_"))
		cat(paste0('         * write file: `lambda_estimate_ecr.txt`'),'\n')
		mu <- read.table("muValues.txt")[max(burn) + Kindex,] # auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
		mu.mcmc <- array(data = NA, dim = c(m,K,p))
		for(k in 1:K){
			mu.mcmc[,k,] <- as.matrix(mu[,k + K*((1:p)-1)])
		}
		mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations$ECR)$output
		write.table(mu.mcmc, file = 'reordered_mu_ecr.txt')
		cat(paste0('         * write file: `reordered_mu_ecr.txt`'),'\n')
		mu.mean <- array(data = NA, dim = c(K,p))
		mu.map <- array(data = NA, dim = c(K,p))
		for(k in 1:K){
			for(j in 1:p){
				mu.mean[k,j] <- mean(mu.mcmc[,k,j])
				mu.map[k,j] <- mu.mcmc[mapIndex,k,j]
			}
		}
		write.table(mu.mean, file = 'mu_estimate_ecr.txt')
		cat(paste0('         * write file: `mu_estimate_ecr.txt`'),'\n')

                SigmaINV <- as.matrix(read.table("sigmainvValues.txt")[max(burn) + Kindex,]) # auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
                SigmaINV.mcmc <- array(data = NA, dim = c(m,K,p))
                for(k in 1:K){
                        SigmaINV.mcmc[,k,] <- as.matrix(SigmaINV[,((k-1)*p + 1):(k*p)])
                }
                SigmaINV.mcmc <- permute.mcmc(SigmaINV.mcmc, ls$permutations$ECR)$output
                Sigma.mcmc <- 1/SigmaINV.mcmc
		write.table(Sigma.mcmc, file = 'reordered_sigma_ecr.txt')
		cat(paste0('         * write file: `reordered_sigma_ecr.txt`'),'\n')
		sigma.mean <- array(data = NA, dim = c(K,p))
		sigma.map <- array(data = NA, dim = c(K,p))
		for(k in 1:K){
			for(j in 1:p){
				sigma.mean[k,j] <- mean(Sigma.mcmc[,k,j])
				sigma.map[k,j] <- Sigma.mcmc[mapIndex,k,j]
			}
		}
		write.table(sigma.mean, file = 'sigma_estimate_ecr.txt')
		cat(paste0('         * write file: `sigma_estimate_ecr.txt`'),'\n')


		w.mcmc <- array(as.matrix(read.table("wValues.txt")[max(burn) + Kindex, ]),dim = c(length(Kindex),K,1))
		w.mcmc_raw <- w.mcmc
		w.mcmc <- array(w.mcmc[,1:K,],dim = c(length(Kindex),K,1))
		w.mcmc <- permute.mcmc(w.mcmc, ls$permutations$"ECR")$output
		w.mean <- numeric(K)
		w.map <- numeric(K)
		for(k in 1:K){
			w.mean[k] <- mean(w.mcmc[,k,1])
			w.map[k] <- w.mcmc[mapIndex,k,1]
		}
		write.table(w.mcmc, file = 'reordered_weights_ecr.txt')
		cat(paste0('         * write file: `reordered_weights_ecr.txt`'),'\n')
		write.table(w.mean, file = 'weights_estimate_ecr.txt')
		cat(paste0('         * write file: `weights_estimate_ecr.txt`'),'\n')
		write.table(ls$clusters, file = "singleBestClusterings.txt", quote = FALSE, row.names = FALSE)
		cat(paste0('         * write file: `singleBestClusterings.txt`'),'\n')
		zMAP <- read.table("singleBestClusterings.txt",header=TRUE)[,1]


		mapAllocationsPosteriorProbs <- numeric(n)
		for(i in 1:n){
			mapAllocationsPosteriorProbs[i] <- length(which(allocationsECR[,i] == zMAP[i]))
		}
		mapAllocationsPosteriorProbs <- mapAllocationsPosteriorProbs/dim(allocationsECR)[1]
		write.table(mapAllocationsPosteriorProbs, file = "classificationProbabilities.txt")
		cat(paste0('         * write file: `classificationProbabilities.txt`'),'\n')
		cat(paste0('-    Done.'), '\n')
	}
	setwd("../")

}






