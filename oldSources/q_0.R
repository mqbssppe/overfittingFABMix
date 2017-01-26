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
                diag(A[k,,]) <- 1/( suff_statistics$cluster_size[k]*diag(SigmaINV[k,,]) + diag(T_INV))
                B[k,] <- SigmaINV[k,,] %*% (suff_statistics$sx[k,] ) + priorConst1
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

### compute only the matrices
compute_A_B_G_D <- function(SigmaINV, suff_statistics, OmegaINV,K, T_INV){
        A <- array(data = 0, dim = c(K,p,p))
        B <- mu <- array(data = 0, dim = c(K,p))
        G <- 0
        D <- 0
        for(k in 1:K){
                diag(A[k,,]) <- 1/(suff_statistics$cluster_size[k]*diag(SigmaINV[k,,]) + diag(T_INV))
        }
        results <- vector("list", length=4)
        results[[1]] <- A
        results[[2]] <- B
        results[[3]] <- G
        results[[4]] <- D
        names(results) <- c("A","B","G","D")
        return(results)
}


# dirichlet function
myDirichlet <- function(alpha){
        k <- length(alpha)
        theta <- rgamma(k, shape = alpha, rate = 1)
        return(theta/sum(theta))
}

#simulate z and mixture weights from standard Gibbs
update_z <- function(w, mu, Lambda, y, SigmaINV, K){
        probs <- array(data = 0, dim =c(n,K))
        for(k in 1:K){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p , byrow = TRUE)
                probs[,k] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% SigmaINV[k,,] %*% tmp) )}) + 0.5*sum(log(diag(SigmaINV[k,,]))) 
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
                x_var <- SigmaINV[k,,]
                probs[index] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
        }
        return(sum(probs))
}


# update SigmaINV
update_SigmaINV <- function(z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
        SigmaINV <- array(data = 0, dim = c(K,p,p))
        for(k in 1:K){
                index <- which(z == k)
                n_k <- length(index)
                s <- numeric(p)
                #s <- array(0, dim = c(p,p))
                for(i in index){
                        tmp <- x_data[i, ] - mu[k,]
                        # apo ton pinaka tmp x tmp xreiazomai mono ta diagonia stoixeia
                        #s <- s + tmp %*% t(tmp)
                        s <- s + tmp^2
                }
                diag(SigmaINV[k,,]) <- rgamma(p,shape = alpha_sigma + n_k/2, rate = beta_sigma + s/2)
        }
        return(SigmaINV)
}




update_SigmaINV_faster <- function(z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
        SigmaINV <- array(data = 0, dim = c(K,p,p))
	alive <- as.numeric(names(table(z)))
	for (k in 1:K){
	        s <- numeric(p)
		ind <- which(z == k)
		n_k <- length(ind)
		if(n_k > 0){
			tmp <- matrix(mu[k, ], nrow = n_k, ncol = p, byrow = TRUE) 
			s <- colSums((x_data[ind, ] - tmp)^2)
		}
	        diag(SigmaINV[k, , ]) <- rgamma(p,shape = alpha_sigma + n_k/2, rate = beta_sigma + s/2)
	}
        return(SigmaINV)
}


#update OmegaINV        #den xreiazetai
update_OmegaINV <- function(Lambda, K, g, h){
        OmegaINV <- array(data = 0, dim = c(q,q))
        betaVector <- numeric(q)
        for(k in 1:K){
                betaVector <- betaVector + colSums(array(Lambda[k,,]^2,dim = c(p,q)))
        }
        diag(OmegaINV) <- rgamma(q,shape = g + K*p/2, rate = h + betaVector/2)
        return(OmegaINV)
}

mh_move_on_z <- function(K, z, y, suff_statistics, SigmaINV, OmegaINV, current_matrix, T_INV, alpha_prior){
	j <- sample(K, 2, replace = FALSE)
	index1 <- which(z == j[1])
	accepted = 0
	if( suff_statistics$cluster_size[j[1]] > 0 ){
		sxx <- array(data = 0, dim = c(p,p))
		syy <- array(data = 0, dim = c(q,q))
		sxy <- array(data = 0, dim = c(p,q))
		n_selected_obs <- floor(runif(1)*suff_statistics$cluster_size[j[1]]) + 1
		selected_obs <- index1[sample(suff_statistics$cluster_size[j[1]], n_selected_obs)]
		proposed_suff_statistics <- suff_statistics
		proposed_suff_statistics$cluster_size[j] <- suff_statistics$cluster_size[j] + c( -n_selected_obs, n_selected_obs)
		x_sub <- array(x_data[selected_obs, ], dim = c(n_selected_obs, p))
		y_sub <- array(y[selected_obs, ], dim = c(n_selected_obs, q))
		sx <- colSums(x_sub)
		sy <- colSums(y_sub)
		for( i in selected_obs){
			sxx <- sxx + x_data[i,] %*% t(x_data[i,])
			syy <- syy + y[i,] %*% t(y[i,])
			sxy <- sxy + x_data[i,] %*% t(y[i,])
		}
		proposed_suff_statistics$sx[j[1],] <- suff_statistics$sx[j[1],] - sx
		proposed_suff_statistics$sx[j[2],] <- suff_statistics$sx[j[2],] + sx

		proposed_suff_statistics$sy[j[1],] <- suff_statistics$sy[j[1],] - sy
		proposed_suff_statistics$sy[j[2],] <- suff_statistics$sy[j[2],] + sy

		proposed_suff_statistics$sxx[j[1],,] <- suff_statistics$sxx[j[1],,] - sxx
		proposed_suff_statistics$sxx[j[2],,] <- suff_statistics$sxx[j[2],,] + sxx

		proposed_suff_statistics$syy[j[1],,] <- suff_statistics$syy[j[1],,] - syy
		proposed_suff_statistics$syy[j[2],,] <- suff_statistics$syy[j[2],,] + syy

		proposed_suff_statistics$sxy[j[1],,] <- suff_statistics$sxy[j[1],,] - sxy
		proposed_suff_statistics$sxy[j[2],,] <- suff_statistics$sxy[j[2],,] + sxy
		# compute new stuff
		proposed_matrix <- compute_A_B_G_D(SigmaINV = SigmaINV, suff_statistics = proposed_suff_statistics, 
					OmegaINV = OmegaINV,K = K, T_INV = T_INV)


		log_ar <- log(suff_statistics$cluster_size[j[1]]) - 
				log(suff_statistics$cluster_size[j[2]] + n_selected_obs) +
					  lgamma(suff_statistics$cluster_size[j[1]] + 1) + 
				          lgamma(suff_statistics$cluster_size[j[2]] + 1) - 
				 lgamma(proposed_suff_statistics$cluster_size[j[1]] + 1) - 
				 lgamma(proposed_suff_statistics$cluster_size[j[2]] + 1)
		log_ar <- log_ar - sum( lgamma(alpha_prior[j] + suff_statistics$cluster_size[j]) ) + 
				   sum( lgamma(alpha_prior[j] + proposed_suff_statistics$cluster_size[j]) ) +
				   0.5*(proposed_suff_statistics$cluster_size[j[1]] - suff_statistics$cluster_size[j[1]])*sum(log(diag(SigmaINV[j[1],,]))) + 
				   0.5*(proposed_suff_statistics$cluster_size[j[2]] - suff_statistics$cluster_size[j[2]])*sum(log(diag(SigmaINV[j[2],,])))
		for(k in j){
			tmp <- sum(diag(SigmaINV[k,,]%*%(proposed_suff_statistics$sxx[k,,] - suff_statistics$sxx[k,,])))
			log_ar <- log_ar - 0.5*tmp
		}
		for(r in 1:p){
			subIndex <- 1:v_r[r]
			for(k in j){
				tmpDET <- det(array(proposed_matrix$G[k,r,subIndex,subIndex],dim = c(v_r[r],v_r[r])))
				#if(tmpDET <= 0){log_ar <- -9999999; print("wtf")}else{
					tmp <- (proposed_matrix$A[k,r,r]*(proposed_suff_statistics$sx[k,r]^2) - 
							current_matrix$A[k,r,r]*(suff_statistics$sx[k,r]^2))*(SigmaINV[k,r,r]^2) + 
						sum( diag(
							proposed_matrix$D[k,r,subIndex]%*% t(proposed_matrix$D[k,r,subIndex]) %*% proposed_matrix$G[k,r,subIndex,subIndex] - 
							current_matrix$D[k,r,subIndex] %*% t(current_matrix$D[k,r,subIndex]) %*% current_matrix$G[k,r,subIndex,subIndex] 
							)) + 
						log(tmpDET/det(array(current_matrix$G[k,r,subIndex,subIndex],dim = c(v_r[r],v_r[r]))))
					log_ar <- log_ar + 0.5*tmp
				#}
			}
		}
		log_ar <- log_ar + 0.5*( log(det(proposed_matrix$A[k,,])) - log(det(current_matrix$A[k,,])) )
		#print(log_ar)
		z_new <- z
		z_new[selected_obs] <- rep(j[2],n_selected_obs)
		#par(mfrow = c(1,2))
		#plot(x_data[,c(1,2)],col = z)
		#plot(x_data[,c(1,2)],col = z_new)
		if( log(runif(1)) < log_ar ){
			z <- z_new
			suff_statistics <- proposed_suff_statistics
			current_matrix <- proposed_matrix
			accepted = 1
		}
		
	}
	results <- vector("list",length = 4)
	results[[1]] <- accepted
	results[[2]] <- z
	results[[3]] <- suff_statistics
	results[[4]] <- current_matrix
	names(results) <- c("accepted","z","suff_statistics","current_matrix")
	return(results)
}

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################

truncatedDirichletProcess <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, progressGraphs, start_values, q){
	if(missing(originalX)){originalX <- x_data}
	if(missing(x_data)){stop('x_data not provided.')}
	if(q != 0){stop('q should be equal to zero')}
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
        SigmaINV.values <- array(data = 0, dim = c(K,p,p))
        Lambda.values <- 0

	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	omega <- OmegaINV.constant
	omega <- 1/omega

	if(start_values == FALSE){
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			diag(SigmaINV.values[k,,]) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
			#mu.values[iter,k,] <- mu.true[k,]
#			for(r in 1:p){
#				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]])
#			}
		}
#		for(i in 1:n){
#			y[i,] <- rnorm(q,mean = 0,sd = 1)
#		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
	}else{
#		cat(paste0('reading starting values... '))	
#		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt')[1,])
#		omega <- OmegaINV.constant
		tmp1 <- read.table("muValues.txt")
		tmp2 <- as.matrix(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			diag(SigmaINV.values[k, , ]) <- as.matrix(tmp2[,((k-1)*p + 1):(k*p)])
#			Lambda.values[k, , ] <- matrix(as.matrix(read.table(paste0("LambdaValues",k,".txt"))),nrow = p, ncol = q, byrow=TRUE) 
		}
#		y <- matrix(as.matrix(read.table('yValues.txt')), nrow = n , ncol = q)
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
#	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
#	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	print(OmegaINV.constant)
	for (iter in 2:m){
#		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
		#print(OmegaINV.constant)
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda(SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV)

		mu.values <- f2$mu

		f2 <- update_z(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
                                SigmaINV =array(SigmaINV.values,dim = c(K,p,p)), K = K)
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
	#			cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
	#			cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					myInd <- which(z == k)
					lMyInd <- length(myInd)
	#				for(r in 1:p){
	#					cat(Lambda.values[k, r, ], " ", file = LambdaConnection[[k]], append = TRUE)
	#				}
					cat(diag(SigmaINV.values[k,,]), " ", file = sigmainvConnection, append = TRUE)
	#				cat('\n', file = LambdaConnection[[k]], append = TRUE)
				}
				cat('\n', file = sigmainvConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	#close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	#close(omegainvConnection)
	close(logLConnection)
#	for(k in 1:K){
#		close(regulonExpressionConnection[[k]])
#		close(LambdaConnection[[k]])
#	}
	setwd("../")

}

library('doParallel')


log_dirichlet_pdf <- function(alpha, weights){
	normConstant <- sum( lgamma(alpha) ) - lgamma( sum(alpha) )
	pdf <- sum( (alpha - 1)*log(weights) ) - normConstant
	return(pdf)
}

heated_chains <- function(dirPriorAlphas, rawData, outDir, Kmax, mCycles, burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize, thinning, nIterPerCycle){
	if(missing(Kmax)){Kmax <- 20}
	if(missing(nIterPerCycle)){nIterPerCycle = 10}
	if( missing(dirPriorAlphas) ){
		nChains <- 8
		dN <- 1
		dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax

#		dirPriorAlphas <- seq(1,10,length = 8)/Kmax
	#	dirPriorAlphas <- c(seq(1, 0.25, length = 4),c((1/Kmax)*seq(3,1.05, length = 4)), c((1/Kmax)*seq(1.04,1, length = 2)))
	#	dirPriorAlphas <- dirPriorAlphas[length(dirPriorAlphas):1]
	}
	nChains <- length(dirPriorAlphas)
	if( range(diff(order(dirPriorAlphas)))[1] != 1){stop('dirPriorAlphas should be in increasing order.')}
	if( range(diff(order(dirPriorAlphas)))[2] != 1){stop('dirPriorAlphas should be in increasing order.')}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 1}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}

	dir.create(outDir)
	setwd(outDir)
	registerDoParallel(cores = nChains)
	outputDirs <- paste0('alpha_',1:nChains)
	originalX <- rawData
	x_data <- originalX
	if( missing(thinning) ){thinning = 1}
	if( thinning < 1 ){ stop('thinning should be larger than or equal to 1.') }
	thinning <- floor(thinning)
	if( missing(normalize) ){normalize <- TRUE}
	cat('\n')
	cat(paste0("-    p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
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
	cat(paste('-    (1) Initializing from priors that lead to overfitting... '))
	d_per_cluster = 4*p 
	initialAlphas <- seq(d_per_cluster/2, d_per_cluster, length = nChains)
	foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		truncatedDirichletProcess(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
			Kmax = Kmax, m = 100, thinning = 1, burn = 99, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
			alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, progressGraphs = FALSE, start_values = FALSE)
	}
	cat(paste(' OK'),'\n')
	cat(paste('-    (2) Initializing the actual model from the previously obtained values... '))
	foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		truncatedDirichletProcess(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
			Kmax = Kmax, m = 300, thinning = 1, burn = 299, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
			alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, progressGraphs = FALSE, start_values = TRUE)
	}
	cat(paste(' OK'),'\n')

	for(myChain in 1:nChains){
		kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
	}
	file_names <- list.files(outputDirs[1])	# the name of files are the same for each folder

	zConnection_target <- file(paste0(getwd(),"/zValues.txt"),open = "w")
	sigmainvConnection_target <- file(paste0(getwd(),"/sigmainvValues.txt"),open = "w")
	muConnection_target <- file(paste0(getwd(),"/muValues.txt"),open = "w")
	wConnection_target <- file(paste0(getwd(),"/wValues.txt"),open = "w")
	cllConnection_target <- file(paste0(getwd(),"/cllValues.txt"),open = "w")
	
	#	main loops
	weights <- array(data = NA, dim = c(2, Kmax))
	#par(mfrow = c(1, 2))
	bb <- nIterPerCycle - 1
	for( iteration in 2:mCycles ){
		
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			truncatedDirichletProcess(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax,  m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
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
		if( (iteration %% 10) == 0 ){
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


