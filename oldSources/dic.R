#just to read the data
#set.seed(1)
#source('sim_full.R')
burn <- 1:1

Q <- 10
nClusters <- c()
for(i in 1:Q){
	cat(paste0("*nFactors = ",i),"\n")
	j <- 1
	lj <- read.table(paste0("q_",i,"/kValues.txt"))
	l <- lj
	nClusters <- cbind(nClusters,l[,1])
}

Km <- length(list.files("q_1",pattern="LambdaValues"))
p <- dim(x_data)[2]
n <- dim(x_data)[1]

library('label.switching')


ct <- -(p/2)*log(2*pi)

observed.log.likelihood0 <- function(w, mu, Lambda, Sigma, z){
        probs <- numeric(n)
	if( is.null(z) ){alive <- 1}else{
        alive <- as.numeric(names(table(z)))}
	newW <- numeric(length(w))
	newW[alive] <- w[alive]/sum(w[alive])
	loggedValues <- array(data = NA, dim = c(n, length(alive)))
	colnames(loggedValues) <- alive
        for(k in alive){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
                x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
                diag(x_var) <- diag(x_var) + Sigma[k,]
                x_var <- try(solve(x_var), TRUE)
		loggedValues[ ,as.character(k)] <- log(newW[k]) - 0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) + ct
        }
	lMax <- apply(loggedValues, 1, max)
	if( length(alive) == 1 ){
		logL <- sum(lMax + log( exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ) ) )
	}else{
		logL <- sum(lMax + log( rowSums(exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ))))
	}
        return(	logL )
}


getStuffForDIC <- function(q){
	cat(paste0('************************************       q = ', q))
	setwd(paste0("q_",q))
	z <- as.matrix(read.table("zValues.txt"))
	logl <- read.table("cllValues.txt")
#	tmp <- read.table("kValues.txt",header = TRUE)[,1]
	tmp  <- apply(z,1,function(y){length(table(y))})
	logl <- cbind(tmp, logl)
	logl <- logl[-burn, ]
	z <- z[-burn,]
	print(table(logl[,1]))
	K <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
	kSelected <- K
	index <- which(logl[,1] == K)
	Kindex <- index
	logl <- logl[index,2]
	z <- z[index,]
	mapIndex <- which(logl == max(logl))[1]
	zPivot <- as.numeric(z[mapIndex,])

	K <- max(z)
	m <- length(logl)
#	ls <- label.switching(method = c("ECR", "ECR"), zpivot = zPivot, z = z, K = K)
	if(K > 1){
		ls <- label.switching(method = c("ECR"), zpivot = zPivot, z = z, K = K)
		write.table(file = "zMAP_ecr.txt", ls$clusters[1,], quote = FALSE, row.names=F, col.names=F)
	}else{
		ECR <- rep(1,m)
		permutations <- vector('list', length = 1)
		permutations[[1]] <- array(data = ECR, dim = c(m,1))
		names(permutations) <- "ECR"
		ls <- vector('list', length = 1)
		names(ls) <- "permutations"
		ls$permutations <- permutations	
	}
	oldLabels <- 1:K
	index <- Kindex

	
	allocationsECR <- z
	for (i in 1:m){
		myPerm <- order(ls$permutations$"ECR"[i,])
		allocationsECR[i,] <- myPerm[z[i,]]
	}



	#Lambda
	l <- as.matrix(read.table(paste0("LambdaValues",1,".txt")))
	J <- dim(l)[2]
	mcmc <- array(data = NA, dim = c(m,K,J))
	for(k in 1:K){
	#	lMean <- apply(l,2,mean)
	#	lMean.matrix <- matrix(lMean,nrow = p, ncol = q, byrow=TRUE)
		l <- as.matrix(read.table(paste0("LambdaValues",k,".txt")))
		l <- l[-burn,]
		mcmc[,k,] <- l[index,]
	}
	lambda.perm.mcmc <- permute.mcmc(mcmc, ls$permutations$ECR)$output
	lambda.mean <- array(data = NA, dim = c(K,p,q))
	lambda.map <- array(data = NA, dim = c(K,p,q))
	for(k in 1:K){
		lMean <- apply(lambda.perm.mcmc[,k,],2,mean)
		lambda.mean[k,,] <- matrix(lMean,nrow = p, ncol = q, byrow=TRUE)
		lambda.map[k,,] <- matrix(lambda.perm.mcmc[mapIndex,k,],nrow = p, ncol = q, byrow=TRUE)
	}

	mu <- read.table("muValues.txt")[max(burn) + Kindex,] # auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
	mu.mcmc <- array(data = NA, dim = c(m,K,p))
	for(k in 1:K){
		mu.mcmc[,k,] <- as.matrix(mu[,k + Km*((1:p)-1)])
	}
	mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations$ECR)$output
	mu.mean <- array(data = NA, dim = c(K,p))
	mu.map <- array(data = NA, dim = c(K,p))
	for(k in 1:K){
		for(j in 1:p){
			mu.mean[k,j] <- mean(mu.mcmc[,k,j])
			mu.map[k,j] <- mu.mcmc[mapIndex,k,j]
		}
	}
	#
	SigmaINV <- as.matrix(read.table("sigmainvValues.txt")[max(burn) + Kindex,]) # auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
	SigmaINV.mcmc <- array(data = NA, dim = c(m,K,p))
	for(k in 1:K){
		SigmaINV.mcmc[,k,] <- as.matrix(SigmaINV[,((k-1)*p + 1):(k*p)])
	}
	SigmaINV.mcmc <- permute.mcmc(SigmaINV.mcmc, ls$permutations$ECR)$output
	SigmaINV.mean <- array(data = NA, dim = c(K,p))
	SigmaINV.map <- array(data = NA, dim = c(K,p))
	for(k in 1:K){
		for(j in 1:p){
			SigmaINV.mean[k,j] <- mean(SigmaINV.mcmc[,k,j])
			SigmaINV.map[k,j] <- SigmaINV.mcmc[mapIndex,k,j]
		}
	}

	Sigma.mcmc <- 1/SigmaINV.mcmc
	Sigma.mean <- array(data = NA, dim = c(K,p))
	Sigma.map <- array(data = NA, dim = c(K,p))
	for(k in 1:K){
		for(j in 1:p){
			Sigma.mean[k,j] <- min(c(1,mean(Sigma.mcmc[,k,j])))
			Sigma.map[k,j] <- min(c(1,Sigma.mcmc[mapIndex,k,j]))
		}
	}

	#SigmaINV.mean <- as.numeric(apply(SigmaINV,2,mean))

	w.mcmc <- array(as.matrix(read.table("wValues.txt")[max(burn) + Kindex, ]),dim = c(length(Kindex),Km,1))
	w.mcmc <- array(w.mcmc[,1:K,],dim = c(length(Kindex),K,1))
	w.mcmc <- permute.mcmc(w.mcmc, ls$permutations$"ECR")$output
	w.mean <- numeric(K)
	w.map <- numeric(K)
	for(k in 1:K){
		w.mean[k] <- mean(w.mcmc[,k,1])
		w.map[k] <- w.mcmc[mapIndex,k,1]
	}




	logL.theta_hat <- observed.log.likelihood0(w = w.mean, mu = mu.mean, Lambda = lambda.mean, Sigma = Sigma.mean, z = ls$clusters[1,])
	logL.theta_map <- observed.log.likelihood0(w = w.map, mu = mu.map, Lambda = lambda.map, Sigma = Sigma.map, z = ls$clusters[1,])
	cll <- 0
	aic <- 0
	bic <- 0
	maxL <- 1
	i <- 1
	lambda.current <- array(data = NA, dim = c(K,p,q))
	Sigma.current <- Sigma.mcmc[i, , ]
	mu.current <- mu.mcmc[i,,]
	if( K == 1 ){ 	
		Sigma.current <- array(Sigma.current, dim = c(1, p))
		mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
	}
	for(k in 1:K){
		ldraw <- lambda.perm.mcmc[i,k,]
		lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
		for(i1 in 1:p){
			if( Sigma.current[k, i1] > 100 ){ Sigma.current[k, i1] <- 100 }
		}
	}
	obsL <- observed.log.likelihood0(w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = allocationsECR[i,])
	maxL <- obsL
	cll <- cll + obsL
	aic <- aic + obsL
	bic <- bic + obsL
	iterMax <- i
	for(i in 2:m){
		lambda.current <- array(data = NA, dim = c(K,p,q))
		Sigma.current <- Sigma.mcmc[i, , ]
		mu.current <- mu.mcmc[i,,]
		if( K == 1 ){ 	
			Sigma.current <- array(Sigma.current, dim = c(1, p))
			mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
		}
		for(k in 1:K){
			ldraw <- lambda.perm.mcmc[i,k,]
			lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
			for(i1 in 1:p){
				if( Sigma.current[k, i1] > 100 ){ Sigma.current[k, i1] <- 100 }
			}
		}
		obsL <- observed.log.likelihood0(w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = allocationsECR[i,])
		if( obsL > maxL ){
			maxL <- obsL	
			iterMax <- i
		}
		cll <- cll + obsL
		aic <- aic + obsL
		bic <- bic + obsL
	}
	cll <- cll/m
	dic_classic <- -4*cll + 2*logL.theta_hat
	dic_star <- -6*cll + 4*logL.theta_hat
#	dic_classicMAP <- -4*cll + 2*logL.theta_map
#	dic_starMAP <- -6*cll + 4*logL.theta_map
	dic_classicMAP <- -4*cll + 2*maxL
	dic_starMAP <- -6*cll + 4*maxL
	Km <- kSelected
	aic <- -2*aic/m + 2*(Km*( 2*p+p*q - q*(q-1)/2 ) + Km - 1 ) 
	bic <- -2*bic/m + log(n)*(Km*( 2*p+p*q - q*(q-1)/2 ) + Km - 1 )
	aic_MAX <- -2*maxL + 2*(Km*( 2*p+p*q - q*(q-1)/2 ) + Km - 1 ) 
	bic_MAX <- -2*maxL + log(n)*(Km*( 2*p+p*q - q*(q-1)/2 ) + Km - 1 )

	dic <- c(aic, bic, dic_classic, dic_star, dic_classicMAP, dic_starMAP, aic_MAX, bic_MAX)
	setwd("../")
	return(dic)
}




#dics2 <- numeric(Q)
#lV1 <- lV2 <- vector("list",length=10)
#for(q in 1:Q){
#	dd <- getStuffForDIC3(q, dirindex = dirIndex[q])
#	dics2[q] <- dd$dic
#	lV1[[q]] <- dd$lV1
#	lV2[[q]] <- dd$lV2
#}

par(mfrow = c(1,1))
library('RColorBrewer')
myPCH <- c(0:7)
myCol <- brewer.pal(8, "Set2")
dics <- array(data = NA, dim = c(Q, 8))
colnames(dics) <- c('AIC', 'BIC', 'DIC', 'DIC*', 'DIC_map', 'DIC*_map', 'AIC_map', 'BIC_map')
for(q in 1:Q){
	dics[q, ] <- getStuffForDIC(q)
	if(q > 1){
		matplot(dics[1:q, ], type = "b", pch = myPCH, col = 1)
		legend('bottomleft', colnames(dics), pch = myPCH)
	}
}




#pdf(file = "informationCriteria.pdf",width = 18, height = 6)
par(mfrow = c(1,2))
subIndex <- c(5,7,8)
boxplot(nClusters[-burn,],xlab = "number of factors", ylab = "nClusters",main = "Posterior of (K|data,q)")
abline(h = K.true,col = 2,lty = 2)
abline(v = q.true,col = 2,lty = 2)
matplot(dics[,subIndex], type = "b", lwd = 2, pch = myPCH[subIndex], col = myCol[subIndex], ylab = 'criterion', xlab = 'number of factors')
for(i in subIndex){
	j <- order(dics[,i])[1]
	points(j, dics[j,i], pch = 1, col = myCol[i], cex = 3)
}
legend('topright', colnames(dics[,subIndex]), pch = myPCH[subIndex], col = myCol[subIndex], lwd = 2)
#dev.off()



#matplot(dics[, ], ylim = c(range(dics[-1,])), type = "b", pch = myPCH, col = 1, ylab = 'criterion', xlab = 'number of factors')
#for(i in 1:dim(dics)[2]){
#	j <- order(dics[,i])[1]
#	points(j, dics[j,i], pch = 1, col = 2, cex = 2)
#}
#legend('topleft', colnames(dics), pch = myPCH)


if(1 > 2){
	library('mclust')
	mod2 = Mclust(x_data, prior = priorControl(functionName="defaultPrior", shrinkage = 0.0), G = 1:15)
	mod1 = Mclust(originalX, prior = priorControl(functionName="defaultPrior", shrinkage = 0.0), G = 1:15)
	plot(mod1, what = "BIC", ylim = range(mod1$BIC[,-(1:2)], na.rm = TRUE),legendArgs = list(x = "bottomright"))
	plot(mod2, what = "BIC", ylim = range(mod2$BIC[,-(1:2)], na.rm = TRUE),legendArgs = list(x = "bottomright"))
	print(mod1$G)
	print(mod2$G)
}



