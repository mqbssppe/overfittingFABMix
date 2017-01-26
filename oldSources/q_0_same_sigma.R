# computation of log{ f(d_il = 0| ...)/f(d_il = 1| ...) }
library(MASS)
#library("ellipse")

#library(mvtnorm)
# MALLON prepei to delta_i na min ginetai pote to mideniko dianisma
update_delta_y_i <- function(x, mu, SigmaINV, Lambda){
	# 2. gibbs for y_i|...
	centre_x <- x - mu
	Alpha <- t(Lambda) %*% SigmaINV %*% Lambda
	diag(Alpha) <- diag(Alpha) + 1
	Alpha <- solve(Alpha)

	y_mean <- Alpha %*% t(Lambda) %*% SigmaINV %*% centre_x
	y_var <- Alpha
	y_new <- mvrnorm(n = 1, mu = y_mean, Sigma = y_var)

	results <- vector("list", length=1)
	results[[1]] <- y_new
	names(results) <- c("y")
	return(results)
}

update_all_y <- function(x, mu, SigmaINV, Lambda, z){
	# 2. gibbs for y_i|...
	n_k <- table(z)
	alive <- as.numeric(names(n_k))
	y_new <- array(data = 0, dim = c(n, q))
	for(k in alive){
		ind <- which(z == k)
		center_x <- x_data[ind, ] - matrix(mu[k,], nrow = as.numeric(n_k[as.character(k)]), ncol = p, byrow=TRUE)
		Alpha <- t(Lambda[k, , ]) %*% SigmaINV %*% Lambda[k, , ]
		diag(Alpha) <- diag(Alpha) + 1
		Alpha <- solve(Alpha)
		tmp <- Alpha %*% t(Lambda[k, , ]) %*% SigmaINV
		y_mean <- t(apply(center_x, 1, function(tk){return(tmp %*% tk)} ))
		if(q == 1){ y_mean <- t(y_mean) }
		y_new[ind, ] <- t( apply( y_mean, 1, function( tk ){ return( mvrnorm(n = 1, mu = tk, Sigma = Alpha) ) } ) )
	}

	results <- vector("list", length=1)
	results[[1]] <- y_new
	names(results) <- c("y")
	return(results)
}




# compute sufficient statistics given y, z, K (and x_data)
compute_sufficient_statistics <- function(y, z, K){
        cluster_size <- numeric(K)
        sx  <- array(data = 0, dim = c(K,p))
        sy  <- 0
        sxx <- array(data = 0, dim = c(K,p,p))
        syy <- 0
        sxy <- 0
        for(k in 1:K){
                index <- which(z == k)
                cluster_size[k] <- length(index)
                if( cluster_size[k] > 0){
                        sx[k,]  <- colSums(array(x_data[index,],dim = c(cluster_size[k],p)))
                        for( i in index){
                                sxx[k,,] <- sxx[k,,] + x_data[i,] %*% t(x_data[i,])
                        }
                }

        }
        results <- vector("list", length=6)
        names(results) <- c("cluster_size","sx","sy","sxx","syy","sxy")
        results[[1]] <- cluster_size
        results[[2]] <- sx
        results[[3]] <- sy
        results[[4]] <- sxx
        results[[5]] <- syy
        results[[6]] <- sxy
        return(results)
}

# compute A, B, G (\Gamma) and D (\Delta) matrices
# 	given   (1) sufficient statistics
#		(2) SigmaINV matrix
#		(3) \Omega matrix
#		(4) T_INV, ksi (fixed prior parameters, not as arguments)
compute_A_B_G_D_and_simulate_mu_Lambda <- function(SigmaINV, suff_statistics, OmegaINV, K, priorConst1, T_INV){
        A <- array(data = 0, dim = c(K,p,p))
        B <- mu <- array(data = 0, dim = c(K,p))
        G <- 0
        D <- 0
#       I will also simulate    (1) Lambda|sufficient statistics, ...
#                               (2) mu|Lambda, sufficient statistics, ...
        mu <- array(data = 0, dim = c(K,p))
        Lambdas <- 0
        for(k in 1:K){
                diag(A[k,,]) <- 1/( suff_statistics$cluster_size[k]*diag(SigmaINV) + diag(T_INV))
                B[k,] <- SigmaINV %*% (suff_statistics$sx[k,] ) + priorConst1
                # this is for simulating mu_k 
                mu_mean <- A[k,,] %*% B[k,]
                mu[k,] <- mvrnorm(n = 1, mu = mu_mean, Sigma = A[k,,])  
        }
        results <- vector("list", length=6)
        results[[1]] <- A
        results[[2]] <- B
        results[[3]] <- G
        results[[4]] <- D
        results[[5]] <- Lambdas
        results[[6]] <- mu
        names(results) <- c("A","B","G","D","Lambdas","mu")
        return(results)
}


# dirichlet function
myDirichlet <- function(alpha){
        k <- length(alpha)
        theta <- rgamma(k, shape = alpha, rate = 1)
        return(theta/sum(theta))
}


# collapsed over y
update_z2 <- function(w, mu, Lambda, SigmaINV, K){
        probs <- array(data = 0, dim =c(n,K))
        for(k in 1:K){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p , byrow = TRUE)
                probs[,k] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% SigmaINV %*% tmp) )}) + 0.5*sum(log(diag(SigmaINV))) 
        }
        probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
        z <- apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
        results <- vector("list", length=2)
        names(results) <- c("w","z")
        results[[1]] <- w
        results[[2]] <- z
        return(results)
}




complete.log.likelihood <- function(w, mu, Lambda, SigmaINV, z){
        probs <- numeric(n)
        alive <- as.numeric(names(table(z)))
        for(k in alive){
                index <- which(z == k)
                center_x <- x_data[index,] - matrix(mu[k,], nrow = length(index), ncol = p, byrow=TRUE)
                x_var <- SigmaINV
                probs[index] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
        }
        return(sum(probs))
}


update_SigmaINV_faster <- function(z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
        SigmaINV <- array(data = 0, dim = c(p,p))
        s <- numeric(p)
	alive <- as.numeric(names(table(z)))
	for (k in alive){
		ind <- which(z == k)
		n_k <- length(ind)
		tmp <- matrix(mu[k, ], nrow = n_k, ncol = p, byrow = TRUE) 
		s <- s + colSums((x_data[ind, ] - tmp)^2)
	}
        diag(SigmaINV) <- rgamma(p,shape = alpha_sigma + n/2, rate = beta_sigma + s/2)
        return(SigmaINV)
}



#update OmegaINV
update_OmegaINV <- function(Lambda, K, g, h){
	OmegaINV <- array(data = 0, dim = c(q,q))
	betaVector <- numeric(q)
	for(k in 1:K){
		betaVector <- betaVector + colSums(array(Lambda[k,,]^2,dim = c(p,q)))
	}
	diag(OmegaINV) <- rgamma(q,shape = g + K*p/2, rate = h + betaVector/2)
	return(OmegaINV)
}

################################################################################################################
################################################################################################################

truncatedDirichletProcess <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, progressGraphs, start_values, q){
	if(missing(originalX)){originalX <- x_data}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 0.001}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 0.001}
	if(missing(progressGraphs)){progressGraphs <- FALSE}
	if(missing(start_values)){start_values <- FALSE}
	if(progressGraphs == TRUE){library("ellipse")}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- rep(1,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
        OmegaINV.values <- 0
        OmegaINV.constant <- 0
        OmegaINV.constantINITIAL <- OmegaINV.constant
        SigmaINV.values <- array(data = 0, dim = c(p,p))
        Lambda.values <- 0
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	iter <- 1
	omega <- OmegaINV.constant
	omega <- 1/omega
	if(start_values == FALSE){
		diag(SigmaINV.values) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
	}else{
#		cat(paste0('reading starting values... '))	
		tmp1 <- read.table("muValues.txt")
		diag(SigmaINV.values) <- as.numeric(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
		}
		w.values <- as.numeric(read.table("wValues.txt"))
		z <- as.numeric(read.table("zValues.txt"))
#		cat(paste0('done.'),'\n')
	}

	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
#	regulonExpressionConnection <- vector('list',length = K)
#	regulonExpression <- array(data = 0, dim = c(K, p, q))
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	print(OmegaINV.constant)
	for (iter in 2:m){
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K)

		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV)

		mu.values <- f2$mu
		f2 <- update_z2(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), 
						SigmaINV = SigmaINV.values, K = K)
		z <- f2$z

		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){
			index <- which(z == k)
			cluster_size[k] <- length(index)
		}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)

		SigmaINV.values <- update_SigmaINV_faster(z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			#cat(paste0(iter," iterations done."),"\n")
			#print(table(z))
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood(w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(SigmaINV.values), file = sigmainvConnection, '\n', append = TRUE)
			}
		}
	}

	close(zConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(logLConnection)
	setwd("../")

}

library('doParallel')


log_dirichlet_pdf <- function(alpha, weights){
	normConstant <- sum( lgamma(alpha) ) - lgamma( sum(alpha) )
	pdf <- sum( (alpha - 1)*log(weights) ) - normConstant
	return(pdf)
}

heated_chains <- function(dirPriorAlphas, rawData, outDir, Kmax, mCycles, burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize, thinning){
	if(missing(Kmax)){Kmax <- 20}
	if( missing(dirPriorAlphas) ){
#		dirPriorAlphas <- c((1/Kmax)*exp(seq(0,2,length = 10)),0.5,1)
		dirPriorAlphas <- (1/Kmax)*exp( seq(0,4,length = 12) )
		#dirPriorAlphas <- (1/Kmax)*exp(seq(0,2,length = 10))
	#	dirPriorAlphas <- c(seq(1, 0.25, length = 4),c((1/Kmax)*seq(3,1.05, length = 4)), c((1/Kmax)*seq(1.04,1, length = 2)))
	#	dirPriorAlphas <- dirPriorAlphas[length(dirPriorAlphas):1]
	}
	nChains <- length(dirPriorAlphas)
	if( range(diff(order(dirPriorAlphas)))[1] != 1){stop('dirPriorAlphas should be in increasing order.')}
	if( range(diff(order(dirPriorAlphas)))[2] != 1){stop('dirPriorAlphas should be in increasing order.')}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 0.001}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 0.001}

	dir.create(outDir)
	setwd(outDir)
	registerDoParallel(cores = nChains)
	outputDirs <- paste0('alpha_',1:nChains)
	originalX <- rawData
	x_data <- originalX
	if( missing(thinning) ){thinning = 10}
	if( thinning < 1 ){ stop('thinning should be larger than or equal to 1.') }
	thinning <- floor(thinning)
	if( missing(normalize) ){normalize <- TRUE}
	cat('\n')
	cat(paste0("-    p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	cat('-    Same covariance matrix per cluster','\n')
	cat(paste0('-    Using Nchains = ', nChains),'\n')
	cat(paste0('-    Target posterior distribution corresponds to alpha = ', dirPriorAlphas[1]),'\n')
	if( normalize == TRUE ){
		x_data <- scale(originalX, center = TRUE, scale = TRUE)
		cat('-    The sampler uses standardized data.','\n')
	}
	if( normalize == FALSE ){
		x_data <- rawData
		cat('-    The sampler uses raw data.','\n')
	}
	kValues <- array(data = NA, dim = c(mCycles, nChains))
	mh_acceptance_rate <- 0
	dir.create('tmpDir')
	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}

	#	initialization
	iteration <- 1
	foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		truncatedDirichletProcess(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
			Kmax = Kmax, m = 10, thinning = 1, burn = 9, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
			alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, progressGraphs = FALSE, start_values = FALSE)
	}

	for(myChain in 1:nChains){
		kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
	}
	file_names <- list.files(outputDirs[1])	# the name of files are the same for each folder

	#	connections for saving the MCMC output corresponding to the target posterior distribution
	zConnection_target <- file(paste0(getwd(),"/zValues.txt"),open = "w")
	sigmainvConnection_target <- file(paste0(getwd(),"/sigmainvValues.txt"),open = "w")
	muConnection_target <- file(paste0(getwd(),"/muValues.txt"),open = "w")
	wConnection_target <- file(paste0(getwd(),"/wValues.txt"),open = "w")
	cllConnection_target <- file(paste0(getwd(),"/cllValues.txt"),open = "w")
	
	#	main loops
	weights <- array(data = NA, dim = c(2, Kmax))
	#par(mfrow = c(1, 2))
	for( iteration in 2:mCycles ){
		
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			truncatedDirichletProcess(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 3, thinning = 1, burn = 2, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, progressGraphs = FALSE, start_values = TRUE)
			kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
		}


		chains <- sample(nChains - 1, 1)
		chains <- c(chains, chains + 1)
		weights[1, ] <- as.numeric(read.table( paste0(outputDirs[ chains[1] ],'/wValues.txt') ))
		weights[2, ] <- as.numeric(read.table( paste0(outputDirs[ chains[2] ],'/wValues.txt') ))

		mh_denom <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[1, ] ) 
				+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[2, ] )
		mh_nom   <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[1, ] ) 
				+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[2, ] )
		mh_ratio <- mh_nom - mh_denom
		if( log(runif(1)) < mh_ratio ){
#		if( 0 < 1 ){

			# dir1 to tmp
			file.copy( 
					from      = paste0(outputDirs[ chains[1] ], '/', file_names), 
					to        = 'tmpDir',
					overwrite = TRUE
				)
			# dir2 to dir1
			file.copy( 
					from      = paste0(outputDirs[ chains[2] ], '/', file_names), 
					to        = outputDirs[ chains[1] ],
					overwrite = TRUE
				)
			# tmp to dir2
			file.copy( 
					from      = paste0('tmpDir', '/', file_names), 
					to        = outputDirs[ chains[2] ],
					overwrite = TRUE
				)
			
			mh_acceptance_rate <- mh_acceptance_rate + 1
		}

		for(myChain in 1:nChains){
			kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
		}

		z <- as.numeric(read.table('alpha_1/zValues.txt'))
		if( (iteration %% 50) == 0 ){
			par(mfrow = c(1,3))
			matplot(kValues[1:iteration, ], type = "l")
			points(1:iteration, kValues[1:iteration, 1], type = "b", col = 1)
			matplot(t(x_data), type = "l", col = as.numeric(as.factor(z)))
			matplot(t(originalX), type = "l", col = as.numeric(as.factor(z)))
			ar <- round(100*mh_acceptance_rate/iteration, 3)
			cat(paste0('-    iteration: ',iteration,'. Metropolis-Hastings acceptance rate: ', ar), '\n')
		}
		if(iteration %% thinning == 0){
			if(iteration > burnCycles){
				w        <- as.numeric(read.table('alpha_1/wValues.txt')) 
				mu       <- as.numeric(read.table('alpha_1/muValues.txt'))
				sigmainv <- as.numeric(read.table('alpha_1/sigmainvValues.txt')) 
				cll      <- as.numeric(read.table('alpha_1/k.and.logl.Values.txt') )[2]
				cat(z       , file = zConnection_target, '\n', append = TRUE)
				cat(w       , file = wConnection_target, '\n', append = TRUE)
				cat(mu      , file = muConnection_target, '\n', append = TRUE)
				cat(sigmainv, file = sigmainvConnection_target, '\n', append = TRUE)
				cat(cll     , file = cllConnection_target, '\n', append = TRUE)
			}
		}
	}
	stopImplicitCluster()
	close(zConnection_target)
	close(wConnection_target)
	close(muConnection_target)
	close(sigmainvConnection_target)
	close(cllConnection_target)
	keepedSeq <- seq(burnCycles + thinning, mCycles, by = thinning)
	if( burnCycles > 0){
		burnedSeq <- 1:burnCycles
		write.table(file = 'burn_in_period_k.txt', kValues[burnedSeq, ])
	}
	write.table(file = 'kValues.txt', kValues[keepedSeq, ], quote = FALSE, row.names = FALSE, col.names = dirPriorAlphas)
	setwd("../")
	cat('-    DONE.','\n')
}


