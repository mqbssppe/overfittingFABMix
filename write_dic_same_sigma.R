#just to read the data
#set.seed(1)
#source('sim_full.R')
burn <- 1:1
library('RColorBrewer')

# change Km as needed
Km <- 20
p <- dim(x_data)[2]
n <- dim(x_data)[1]

library('label.switching')


ct <- -(p/2)*log(2*pi)

observed.log.likelihood0 <- function(w, mu, Lambda, Sigma, z){
        probs <- numeric(n)
	alive <- 1:Km
	newW <- numeric(length(w))
	newW[alive] <- w[alive]/sum(w[alive])
	loggedValues <- array(data = NA, dim = c(n, length(alive)))
	colnames(loggedValues) <- alive
        for(k in alive){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
                x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
                diag(x_var) <- diag(x_var) + Sigma
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


getStuffForDIC <- function(outputFolder, q){
	cat(paste0('-    (4) Computing information criteria for q = ', q), '\n')
	setwd(outputFolder)
	cat(paste0('         - Entering directory: ', getwd()),'\n')
	z <- as.matrix(read.table("zValues.txt"))
	logl <- read.table("cllValues.txt")
        tmp  <- apply(z,1,function(y){length(table(y))})
        logl <- cbind(tmp, logl)
	z <- z[-burn,]
	cat(paste0('            Nclusters:    ', paste(as.character(names(table(logl[,1]))), collapse="     ") ), '\n')
	cat(paste0('            Frequency:    ', paste(as.character(as.numeric(table(logl[,1]))), collapse="    ") ), '\n')
	K <- Km
	kSelected <- K
	index <- 1:dim(z)[1]
	Kindex <- index
	m <- length(index)

	ECR <- matrix(1:Km, nrow = m, ncol = Km, byrow=T)
	permutations <- vector('list', length = 1)
	permutations[[1]] <- ECR
	names(permutations) <- "ECR"
	ls <- vector('list', length = 1)
	names(ls) <- "permutations"
	ls$permutations <- permutations	



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
	for(k in 1:K){
		lMean <- apply(lambda.perm.mcmc[,k,],2,mean)
	}

	mu <- read.table("muValues.txt")[max(burn) + Kindex,] # auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
	mu.mcmc <- array(data = NA, dim = c(m,K,p))
	for(k in 1:K){
		mu.mcmc[,k,] <- as.matrix(mu[,k + Km*((1:p)-1)])
	}
	mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations$ECR)$output
	#
	SigmaINV <- as.matrix(read.table("sigmainvValues.txt")[max(burn) + Kindex,]) # auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
	SigmaINV.mcmc <- SigmaINV

	Sigma.mcmc <- 1/SigmaINV.mcmc

	#SigmaINV.mean <- as.numeric(apply(SigmaINV,2,mean))

	w.mcmc <- array(as.matrix(read.table("wValues.txt")[max(burn) + Kindex, ]),dim = c(length(Kindex),Km,1))
	w.mcmc <- array(w.mcmc[,1:K,],dim = c(length(Kindex),K,1))
	w.mcmc <- permute.mcmc(w.mcmc, ls$permutations$"ECR")$output


	cll <- 0
	aic <- 0
	bic <- 0
	maxL <- 1
	i <- 1
	lambda.current <- array(data = NA, dim = c(K,p,q))
	Sigma.current <- Sigma.mcmc[i, ]
	mu.current <- mu.mcmc[i,,]
	if( K == 1 ){ 	
		Sigma.current <- array(Sigma.current, dim = c(1, p))
		mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
	}
	for(k in 1:K){
		ldraw <- lambda.perm.mcmc[i,k,]
		lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
		for(i1 in 1:p){
			if( Sigma.current[ i1] > 100 ){ Sigma.current[i1] <- 100 }
		}
	}
	obsL <- observed.log.likelihood0(w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
	maxL <- obsL
	cll <- cll + obsL
	aic <- aic + obsL
	bic <- bic + obsL
	iterMax <- i
	for(i in 2:m){
		lambda.current <- array(data = NA, dim = c(K,p,q))
		Sigma.current <- Sigma.mcmc[i, ]
		mu.current <- mu.mcmc[i,,]
		if( K == 1 ){ 	
			Sigma.current <- array(Sigma.current, dim = c(1, p))
			mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
		}
		for(k in 1:K){
			ldraw <- lambda.perm.mcmc[i,k,]
			lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
			for(i1 in 1:p){
				if( Sigma.current[ i1] > 100 ){ Sigma.current[ i1] <- 100 }
			}
		}
		obsL <- observed.log.likelihood0(w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
		if( obsL > maxL ){
			maxL <- obsL	
			iterMax <- i
		}
		cll <- cll + obsL
		aic <- aic + obsL
		bic <- bic + obsL
	}
	cll <- cll/m
#	dic_classic <- -4*cll + 2*logL.theta_hat
#	dic_star <- -6*cll + 4*logL.theta_hat
#	dic_classicMAP <- -4*cll + 2*logL.theta_map
#	dic_starMAP <- -6*cll + 4*logL.theta_map
	dic_classicMAP <- -4*cll + 2*maxL
	dic_starMAP <- -6*cll + 4*maxL
	dic_classic <- dic_classicMAP
	dic_star <- dic_starMAP
	Km <- kSelected
	aic <- -2*aic/m + 2*(Km*( p+p*q - q*(q-1)/2 ) + p + Km - 1 ) 
	bic <- -2*bic/m + log(n)*(Km*( p+p*q - q*(q-1)/2 ) + p + Km - 1 )
	aic_MAX <- -2*maxL + 2*(Km*( p+p*q - q*(q-1)/2 ) + p + Km - 1 ) 
	bic_MAX <- -2*maxL + log(n)*(Km*( p+p*q - q*(q-1)/2 ) + p + Km - 1 )

	dic <- c(aic, bic, dic_classic, dic_star, dic_classicMAP, dic_starMAP, aic_MAX, bic_MAX)
	names(dic) <- c('AIC', 'BIC', 'DIC1', 'DIC*2', 'DIC', 'DIC_2', 'AIC_map', 'BIC_map')	
	write.table(file = 'informationCriteria.txt', dic[c(1,2,5,6)], col.names = paste0('q_',q), quote = FALSE)
	setwd("../")
	cat(paste0('         - Information criteria written to `', outputFolder,'/informationCriteria.txt`.'), '\n')
}



