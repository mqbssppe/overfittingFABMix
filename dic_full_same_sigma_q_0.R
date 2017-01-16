#just to read the data
#set.seed(1)
#source('sim_full.R')
#burn <- 1:100


# outputDirectory <- 
library('label.switching')

i <- 0
cat(paste0("*nFactors = ",i),"\n")
j <- 1
lj <- read.table(paste0(outputDirectory,"/kValues.txt"), header = TRUE)
l <- lj
nClusters0 <- l[,1]

ct <- -(p/2)*log(2*pi)

observed.log.likelihood0 <- function(w, mu, Sigma, z){
        probs <- numeric(n)
	alive <- 1:Km
	newW <- numeric(length(w))
	newW[alive] <- w[alive]/sum(w[alive])
	loggedValues <- array(data = NA, dim = c(n, length(alive)))
	colnames(loggedValues) <- alive
        for(k in alive){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
                x_var <- array(data = 0, dim = c(p,p))
                diag(x_var) <- Sigma
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


getStuffForDIC <- function(outputDirectory){
##	if(q != 0){stop('q should be zero.')}
#	cat(paste0('************************************       q = ', q))
	setwd(outputDirectory)
	z <- as.matrix(read.table("zValues.txt"))
	logl <- read.table("cllValues.txt")
        tmp  <- apply(z,1,function(y){length(table(y))})
        logl <- cbind(tmp, logl)
	z <- z[-burn,]
	print(table(logl[,1]))
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



	mu <- read.table("muValues.txt")[max(burn) + Kindex,] # auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
	mu.mcmc <- array(data = NA, dim = c(m,K,p))
	for(k in 1:K){
		mu.mcmc[,k,] <- as.matrix(mu[,k + Km*((1:p)-1)])
	}
	mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations$ECR)$output
	#
	SigmaINV <- as.matrix(read.table("sigmainvValues.txt")[max(burn) + Kindex,]) 
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
	Sigma.current <- Sigma.mcmc[i, ]
	mu.current <- mu.mcmc[i,,]
	if( K == 1 ){ 	
		Sigma.current <- array(Sigma.current, dim = c(1, p))
		mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
	}
	for(k in 1:K){
		for(i1 in 1:p){
			if( Sigma.current[ i1] > 100 ){ Sigma.current[ i1] <- 100 }
		}
	}
	obsL <- observed.log.likelihood0(w = w.mcmc[i,,1], mu = mu.current, Sigma = Sigma.current, z = z[i,])
	maxL <- obsL
	cll <- cll + obsL
	aic <- aic + obsL
	bic <- bic + obsL
	iterMax <- i
	for(i in 2:m){
		Sigma.current <- Sigma.mcmc[i, ]
		mu.current <- mu.mcmc[i,,]
		if( K == 1 ){ 	
			Sigma.current <- array(Sigma.current, dim = c(1, p))
			mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
		}
		for(k in 1:K){
			for(i1 in 1:p){
				if( Sigma.current[ i1] > 100 ){ Sigma.current[ i1] <- 100 }
			}
		}
		obsL <- observed.log.likelihood0(w = w.mcmc[i,,1], mu = mu.current, Sigma = Sigma.current, z = z[i,])
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
	aic <- -2*aic/m + 2*(Km*p + p + Km - 1 ) 
	bic <- -2*bic/m + log(n)*(Km*p + p + Km - 1 )
	aic_MAX <- -2*maxL + 2*(Km*p + p + Km - 1 ) 
	bic_MAX <- -2*maxL + log(n)*(Km*p + p + Km - 1 )

	dic <- c(aic, bic, dic_classic, dic_star, dic_classicMAP, dic_starMAP, aic_MAX, bic_MAX, cll)
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

dics0 <- numeric(9)
names(dics0) <- c('AIC', 'BIC', 'DIC', 'DIC*', 'DIC_map', 'DIC*_map', 'AIC_map', 'BIC_map', 'logL')
dics0 <- getStuffForDIC(outputDirectory = outputDirectory)
write.table(dics0, file = paste0(outputDirectory,'/dics0.txt'))


