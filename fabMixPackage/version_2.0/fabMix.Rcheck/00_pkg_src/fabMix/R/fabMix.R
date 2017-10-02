#library('MASS')
#library('mvtnorm')
#library('foreach')
#library('doParallel')
#library('label.switching')

update_all_y <- function(x_data, mu, SigmaINV, Lambda, z){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	n <- dim(x_data)[1]
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

update_all_y_Sj <- function(x_data, mu, SigmaINV, Lambda, z){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	n <- dim(x_data)[1]
	# 2. gibbs for y_i|...
	n_k <- table(z)
	alive <- as.numeric(names(n_k))
	y_new <- array(data = 0, dim = c(n, q))
	for(k in alive){
		ind <- which(z == k)
		center_x <- x_data[ind, ] - matrix(mu[k,], nrow = as.numeric(n_k[as.character(k)]), ncol = p, byrow=TRUE)
		Alpha <- t(Lambda[k, , ]) %*% SigmaINV[k, ,] %*% Lambda[k, , ]
		diag(Alpha) <- diag(Alpha) + 1
		Alpha <- solve(Alpha)
		tmp <- Alpha %*% t(Lambda[k, , ]) %*% SigmaINV[k,,]
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
compute_sufficient_statistics <- function(y, z, K, x_data){
	cluster_size <- numeric(K)
	p <- dim(x_data)[2]
	q <- dim(y)[2]
	sx  <- array(data = 0, dim = c(K,p))
	sy  <- array(data = 0, dim = c(K,q))
	sxx <- array(data = 0, dim = c(K,p,p))
	syy <- array(data = 0, dim = c(K,q,q))
	sxy <- array(data = 0, dim = c(K,p,q))
	for(k in 1:K){
		index <- which(z == k)
		cluster_size[k] <- length(index)
		if( cluster_size[k] > 0){
			sx[k,]  <- colSums(array(x_data[index,],dim = c(cluster_size[k],p)))
			sy[k,]  <- colSums(array(y[index,],dim = c(cluster_size[k],q)))
			for( i in index){
				sxx[k,,] <- sxx[k,,] + x_data[i,] %*% t(x_data[i,])
				syy[k,,] <- syy[k,,] + y[i,] %*% t(y[i,])
				sxy[k,,] <- sxy[k,,] + x_data[i,] %*% t(y[i,])
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
compute_A_B_G_D_and_simulate_mu_Lambda <- function(SigmaINV, suff_statistics, OmegaINV, K, priorConst1, T_INV, v_r){
	p <- dim( suff_statistics$sx )[2]
	q <- dim( suff_statistics$sy )[2]
	A <- array(data = 0, dim = c(K,p,p))
	B <- mu <- array(data = 0, dim = c(K,p))
	G <- array(data = 0, dim = c(K,p,q,q))
	D <- array(data = 0, dim = c(K,p,q))
#	I will also simulate    (1) Lambda|sufficient statistics, ...
#				(2) mu|Lambda, sufficient statistics, ...
	mu <- array(data = 0, dim = c(K,p))
	Lambdas <- array(data = 0, dim = c(K,p,q))
	for(k in 1:K){
		diag(A[k,,]) <- 1/(suff_statistics$cluster_size[k]*diag(SigmaINV) + diag(T_INV))
		for (r in 1:p){
			# computing G[k,r,,]: the actual dimension is v_r x v_r, r = 1,...,p
			subIndex <- 1:v_r[r]
			tmp <- OmegaINV[subIndex,subIndex] + suff_statistics$syy[k,subIndex,subIndex]*SigmaINV[r,r] - 
				A[k,r,r]*(SigmaINV[r,r]^2)*suff_statistics$sy[k,subIndex]%*%t(suff_statistics$sy[k,subIndex])
			G[k,r,subIndex,subIndex] <- solve(tmp)
			#computing D[k,r,]: the actual dimension is v_r, r = 1,...,p
			D[k,r,subIndex] <- suff_statistics$sxy[k,r,subIndex]*SigmaINV[r,r] - A[k,r,r]*(SigmaINV[r,r]^2)*suff_statistics$sx[k,r]*suff_statistics$sy[k,subIndex]
			# this is for simulating Lambdas[k,r,,]
			lambda_mean <- G[k,r,subIndex,subIndex] %*% D[k,r,subIndex]
			Lambdas[k,r,subIndex] <- mvrnorm(n = 1, mu = lambda_mean, Sigma = G[k,r,subIndex,subIndex])
		}
		B[k,] <- SigmaINV %*% (suff_statistics$sx[k,] - array(Lambdas[k,,],dim=c(p,q)) %*% array(suff_statistics$sy[k,], dim=c(q,1)) ) + priorConst1
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

compute_A_B_G_D_and_simulate_mu_Lambda_Sj <- function(SigmaINV, suff_statistics, OmegaINV, K, priorConst1, T_INV, v_r){
	p <- dim( suff_statistics$sx )[2]
	q <- dim( suff_statistics$sy )[2]
	A <- array(data = 0, dim = c(K,p,p))
	B <- mu <- array(data = 0, dim = c(K,p))
	G <- array(data = 0, dim = c(K,p,q,q))
	D <- array(data = 0, dim = c(K,p,q))
#	I will also simulate    (1) Lambda|sufficient statistics, ...
#				(2) mu|Lambda, sufficient statistics, ...
	mu <- array(data = 0, dim = c(K,p))
	Lambdas <- array(data = 0, dim = c(K,p,q))
	for(k in 1:K){
		diag(A[k,,]) <- 1/(suff_statistics$cluster_size[k]*diag(SigmaINV[k,,]) + diag(T_INV))
		for (r in 1:p){
			# computing G[k,r,,]: the actual dimension is v_r x v_r, r = 1,...,p
			subIndex <- 1:v_r[r]
			tmp <- OmegaINV[subIndex,subIndex] + suff_statistics$syy[k,subIndex,subIndex]*SigmaINV[k,r,r] - 
				A[k,r,r]*(SigmaINV[k,r,r]^2)*suff_statistics$sy[k,subIndex]%*%t(suff_statistics$sy[k,subIndex])
			G[k,r,subIndex,subIndex] <- solve(tmp)
			#computing D[k,r,]: the actual dimension is v_r, r = 1,...,p
			D[k,r,subIndex] <- suff_statistics$sxy[k,r,subIndex]*SigmaINV[k,r,r] - A[k,r,r]*(SigmaINV[k,r,r]^2)*suff_statistics$sx[k,r]*suff_statistics$sy[k,subIndex]
			# this is for simulating Lambdas[k,r,,]
			lambda_mean <- G[k,r,subIndex,subIndex] %*% D[k,r,subIndex]
			Lambdas[k,r,subIndex] <- mvrnorm(n = 1, mu = lambda_mean, Sigma = G[k,r,subIndex,subIndex])
		}
		B[k,] <- SigmaINV[k,,] %*% (suff_statistics$sx[k,] - array(Lambdas[k,,],dim=c(p,q)) %*% array(suff_statistics$sy[k,], dim=c(q,1)) ) + priorConst1
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


#simulate z and mixture weights from standard Gibbs
update_z_b <- function(w, mu, Lambda, y, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	q <- dim(y)[2]
	probs <- array(data = 0, dim =c(n,K))
	x_var <- array(data = 0, dim = c(p,p))
	diag(x_var) <- 1/diag(SigmaINV)
	for(k in 1:K){
		center_x <- x_data - t(apply(y,1,function(tmp){return(mu[k,] + array(Lambda[k,,], dim = c(p,q)) %*% tmp)}))
		probs[,k] <- log(w[k])  + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
	}
	probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
	z <- apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
	results <- vector("list", length=2)
	names(results) <- c("w","z")
	results[[1]] <- w
	results[[2]] <- z
	return(results)
}

update_z_b_Sj <- function(w, mu, Lambda, y, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	q <- dim(y)[2]
	probs <- array(data = 0, dim =c(n,K))
	for(k in 1:K){
		x_var <- array(data = 0, dim = c(p,p))
		diag(x_var) <- 1/diag(SigmaINV[k,,])
		center_x <- x_data - t(apply(y,1,function(tmp){return(mu[k,] + array(Lambda[k,,], dim = c(p,q)) %*% tmp)}))
		probs[,k] <- log(w[k])  + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
	}
	probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
	z <- apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
	results <- vector("list", length=2)
	names(results) <- c("w","z")
	results[[1]] <- w
	results[[2]] <- z
	return(results)
}


# using dmvnorm from package mvtnorm
update_z4 <- function(w, mu, Lambda, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	probs <- array(data = 0, dim =c(n,K))
	for(k in 1:K){
		center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
		x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
		diag(x_var) <- diag(x_var) + 1/diag(SigmaINV)
		lpdf <- log(w[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
		probs[,k] <- lpdf
	}
	probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
	z <- apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
	results <- vector("list", length=2)
	names(results) <- c("w","z")
	results[[1]] <- w
	results[[2]] <- z
	return(results)
}

update_z4_Sj <- function(w, mu, Lambda, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	probs <- array(data = 0, dim =c(n,K))
	for(k in 1:K){
		center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
		x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
		diag(x_var) <- diag(x_var) + 1/diag(SigmaINV[k,,])
		lpdf <- log(w[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
		probs[,k] <- lpdf
	}
	probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
	z <- apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
	results <- vector("list", length=2)
	names(results) <- c("w","z")
	results[[1]] <- w
	results[[2]] <- z
	return(results)
}



# using the matrix inversion lemma
update_z2 <- function(w, mu, Lambda, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
        probs <- array(data = 0, dim =c(n,K))
        for(k in 1:K){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
		x_var <- t(Lambda[k,,]) %*% SigmaINV %*% Lambda[k,,]
		diag(x_var) <- diag(x_var) + 1
		x_var <- try(solve(x_var), TRUE)
                if(is.numeric(x_var) == TRUE){
			x_var <- Lambda[k,,] %*% x_var %*% t(Lambda[k,,])
			x_var <- SigmaINV %*% x_var %*% SigmaINV
			x_var <- SigmaINV - x_var
                        probs[,k] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
                }
        }
        probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
        z <- apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
        results <- vector("list", length=2)
        names(results) <- c("w","z")
        results[[1]] <- w
        results[[2]] <- z
        return(results)
}


update_z2_Sj <- function(w, mu, Lambda, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
        probs <- array(data = 0, dim =c(n,K))
        for(k in 1:K){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
		x_var <- t(Lambda[k,,]) %*% SigmaINV[k,,] %*% Lambda[k,,]
		diag(x_var) <- diag(x_var) + 1
		x_var <- try(solve(x_var), TRUE)
                if(is.numeric(x_var) == TRUE){
			x_var <- Lambda[k,,] %*% x_var %*% t(Lambda[k,,])
			x_var <- SigmaINV[k,,] %*% x_var %*% SigmaINV[k,,]
			x_var <- SigmaINV[k,,] - x_var
                        probs[,k] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
                }
        }
        probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
        z <- apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
        results <- vector("list", length=2)
        names(results) <- c("w","z")
        results[[1]] <- w
        results[[2]] <- z
        return(results)
}





complete.log.likelihood <- function(x_data, w, mu, Lambda, SigmaINV, z){
	n <- length(z)
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]

	probs <- numeric(n)
	alive <- as.numeric(names(table(z)))
	for(k in alive){
		index <- which(z == k)
		center_x <- x_data[index,] - matrix(mu[k,], nrow = length(index), ncol = p, byrow=TRUE)
		x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
		diag(x_var) <- diag(x_var) + 1/diag(SigmaINV)
		#x_var <- solve(x_var)
		x_var <- try(solve(x_var), TRUE)
		if(is.numeric(x_var) == TRUE){
			probs[index] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
		}
	}
	return(sum(probs))
}


complete.log.likelihood_Sj <- function(x_data, w, mu, Lambda, SigmaINV, z){
	n <- length(z)
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]

	probs <- numeric(n)
	alive <- as.numeric(names(table(z)))
	for(k in alive){
		index <- which(z == k)
		center_x <- x_data[index,] - matrix(mu[k,], nrow = length(index), ncol = p, byrow=TRUE)
		x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
		diag(x_var) <- diag(x_var) + 1/diag(SigmaINV[k,,])
		#x_var <- solve(x_var)
		x_var <- try(solve(x_var), TRUE)
		if(is.numeric(x_var) == TRUE){
			probs[index] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
		}
	}
	return(sum(probs))
}



update_SigmaINV_faster <- function(x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	n <- length(z)

        SigmaINV <- array(data = 0, dim = c(p,p))
        s <- numeric(p)
	alive <- as.numeric(names(table(z)))
	for (k in alive){
		ind <- which(z == k)
		n_k <- length(ind)
		tmp <- matrix(mu[k, ], nrow = n_k, ncol = p, byrow = TRUE) + t(
			apply
			(
				array(y[ind, ], dim = c(n_k,q)), 1, 
					function(tk)
					{ 
						return(array(Lambda[ k, , ], dim = c(p,q)) %*% tk) 
					}
			)
		)
		s <- s + colSums((x_data[ind, ] - tmp)^2)
	}
        diag(SigmaINV) <- rgamma(p,shape = alpha_sigma + n/2, rate = beta_sigma + s/2)
        return(SigmaINV)
}



update_SigmaINV_faster_Sj <- function(x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	n <- length(z)

        SigmaINV <- array(data = 0, dim = c(K,p,p))
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




#update OmegaINV
update_OmegaINV <- function(Lambda, K, g, h){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
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

overfittingMFA <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
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
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
#	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- diag(var(x_data))
#	diag(T_INV) <- rep(1,p)
#	diag(T_INV) <- rep(100,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		diag(SigmaINV.values) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
			}
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt')[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt")
		diag(SigmaINV.values) <- as.numeric(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			Lambda.values[k, , ] <- matrix(as.matrix(read.table(paste0("LambdaValues",k,".txt"))),nrow = p, ncol = q, byrow=TRUE) 
		}
		y <- matrix(as.matrix(read.table('yValues.txt')), nrow = n , ncol = q)
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
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
#	regulonExpressionConnection <- vector('list',length = K)
	LambdaConnection <- vector('list',length = K)
	for(k in 1:K){
#		regulonExpressionConnection[[k]] <- file(paste0("regulonExpressionValues",k,".txt"),open = "w")  #K x p x q per iteration
		LambdaConnection[[k]] <- file(paste0("LambdaValues",k,".txt"),open = "w")  #K x p x q per iteration
	}
#	regulonExpression <- array(data = 0, dim = c(K, p, q))
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	u_v <- runif(1)
	for (iter in 2:m){
		
#		1
		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K, x_data = x_data)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}else{
			f2 <- update_z_b(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y(x_data = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster(x_data = x_data, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				cat(diag(SigmaINV.values), file = sigmainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					myInd <- which(z == k)
					lMyInd <- length(myInd)
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection[[k]], append = TRUE)
					}
					cat('\n', file = LambdaConnection[[k]], append = TRUE)
				}
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	for(k in 1:K){
#		close(regulonExpressionConnection[[k]])
		close(LambdaConnection[[k]])
	}
	setwd("../")

}




overfittingMFA_Sj <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
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
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
#	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- diag(var(x_data))
#	diag(T_INV) <- rep(1,p)
#	diag(T_INV) <- rep(100,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(K,p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
			}
			diag(SigmaINV.values[k,,]) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt')[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt")
		tmp2 <- as.matrix(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			Lambda.values[k, , ] <- matrix(as.matrix(read.table(paste0("LambdaValues",k,".txt"))),nrow = p, ncol = q, byrow=TRUE)
			diag(SigmaINV.values[k, , ]) <- as.matrix(tmp2[,((k-1)*p + 1):(k*p)]) 
		}
		y <- matrix(as.matrix(read.table('yValues.txt')), nrow = n , ncol = q)
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
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
#	regulonExpressionConnection <- vector('list',length = K)
	LambdaConnection <- vector('list',length = K)
	for(k in 1:K){
#		regulonExpressionConnection[[k]] <- file(paste0("regulonExpressionValues",k,".txt"),open = "w")  #K x p x q per iteration
		LambdaConnection[[k]] <- file(paste0("LambdaValues",k,".txt"),open = "w")  #K x p x q per iteration
	}
#	regulonExpression <- array(data = 0, dim = c(K, p, q))
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	u_v <- runif(1)
	for (iter in 2:m){
		
#		1
		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K, x_data = x_data)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_Sj(SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}else{
			f2 <- update_z_b_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), K = K, x_data = x_data)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y_Sj(x_data = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster_Sj(x_data = x_data, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood_Sj(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					myInd <- which(z == k)
					lMyInd <- length(myInd)
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection[[k]], append = TRUE)
					}
					cat(diag(SigmaINV.values[k,,]), " ", file = sigmainvConnection, append = TRUE)
					cat('\n', file = LambdaConnection[[k]], append = TRUE)
				}
				cat('\n', file = sigmainvConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	for(k in 1:K){
#		close(regulonExpressionConnection[[k]])
		close(LambdaConnection[[k]])
	}
	setwd("../")

}




overfittingMFA_missing_values <- function(missing_entries, x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z){
	# missing_entries: list which contains the row number (1st entry) and column indexes (subsequent entries) for every row containing missing values

	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
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
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data, na.rm = TRUE))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colMeans(x_data, na.rm = TRUE)
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		diag(SigmaINV.values) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
			}
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
		# start from colmeans
		myColMeans <- colMeans(x_data, na.rm = TRUE)
		x_complete <- x_data
		xReplacedValues <- lapply(missing_entries, function(y){ 
					x_complete[y[1], y[-1]] <<- myColMeans[y[-1]] 
				}
			)
		x_mean <- x_complete
		
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt')[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt")
		diag(SigmaINV.values) <- as.numeric(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			Lambda.values[k, , ] <- matrix(as.matrix(read.table(paste0("LambdaValues",k,".txt"))),nrow = p, ncol = q, byrow=TRUE) 
		}
		y <- matrix(as.matrix(read.table('yValues.txt')), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt"))
		z <- as.numeric(read.table("zValues.txt"))
#		cat(paste0('done.'),'\n')
		x_complete <- as.matrix(read.table("x_complete.txt", header=TRUE))
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
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	xConnection <- file("x_complete.txt",open = "w")
	LambdaConnection <- vector('list',length = K)
	for(k in 1:K){
		LambdaConnection[[k]] <- file(paste0("LambdaValues",k,".txt"),open = "w")  #K x p x q per iteration
	}
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
	for (iter in 2:m){
#		1
		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K, x_data = x_complete)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_complete)
		}else{
			f2 <- update_z_b(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = SigmaINV.values, K = K, x_data = x_complete)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y(x_data = x_complete, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster(x_data = x_complete, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)
#		6 {missing data}
		mySD <- 1/sqrt(diag(SigmaINV.values))
		xReplacedValues <- lapply(missing_entries, function(f){ 
					nMissing <- length(f) - 1
					j <- z[f[1]]
					myMean <- mu.values[j, f[-1]] + Lambda.values[j,f[-1], ] %*% t(t(y[f[1], ]))
					x_complete[f[1], f[-1]] <<- rnorm(n = nMissing, mean = myMean, sd = mySD[f[-1]]) 
				}
			)
		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				cat(diag(SigmaINV.values), file = sigmainvConnection, '\n', append = TRUE)
				write.table(x_complete, file = xConnection, row.names = FALSE, quote=F, append=TRUE)
				for(k in 1:K){
					myInd <- which(z == k)
					lMyInd <- length(myInd)
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection[[k]], append = TRUE)
					}
					cat('\n', file = LambdaConnection[[k]], append = TRUE)
				}
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(xConnection)
	for(k in 1:K){
		close(LambdaConnection[[k]])
	}
	setwd("../")

}



overfittingMFA_Sj_missing_values <- function(missing_entries, x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
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
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data, na.rm = TRUE))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colMeans(x_data, na.rm = TRUE)
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(K,p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
			}
			diag(SigmaINV.values[k,,]) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
		# start from colmeans
		myColMeans <- colMeans(x_data, na.rm = TRUE)
		x_complete <- x_data
		xReplacedValues <- lapply(missing_entries, function(y){ 
					x_complete[y[1], y[-1]] <<- myColMeans[y[-1]] 
				}
			)
		x_mean <- x_complete
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt')[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt")
		tmp2 <- as.matrix(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			Lambda.values[k, , ] <- matrix(as.matrix(read.table(paste0("LambdaValues",k,".txt"))),nrow = p, ncol = q, byrow=TRUE)
			diag(SigmaINV.values[k, , ]) <- as.matrix(tmp2[,((k-1)*p + 1):(k*p)]) 
		}
		y <- matrix(as.matrix(read.table('yValues.txt')), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt"))
		z <- as.numeric(read.table("zValues.txt"))
#		cat(paste0('done.'),'\n')
		x_complete <- as.matrix(read.table("x_complete.txt", header=TRUE))
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
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	xConnection <- file("x_complete.txt",open = "w")
	LambdaConnection <- vector('list',length = K)
	for(k in 1:K){
		LambdaConnection[[k]] <- file(paste0("LambdaValues",k,".txt"),open = "w")  #K x p x q per iteration
	}
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
	mySD <- array(data = 0, dim = c(K,p))
	for (iter in 2:m){
	
#		1
		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K, x_data = x_complete)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_Sj(SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_complete)
		}else{
			f2 <- update_z_b_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), K = K, x_data = x_complete)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y_Sj(x_data = x_complete, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster_Sj(x_data = x_complete, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)
#		6 {missing data}
		for(k in 1:K){
			mySD[k,] <-  1/sqrt(diag(SigmaINV.values[k, , ]))
		}
		xReplacedValues <- lapply(missing_entries, function(f){ 
					nMissing <- length(f) - 1
					j <- z[f[1]]
					myMean <- mu.values[j, f[-1]] + Lambda.values[j,f[-1], ] %*% t(t(y[f[1], ]))
					x_complete[f[1], f[-1]] <<- rnorm(n = nMissing, mean = myMean, sd = mySD[j, f[-1]]) 
				}
			)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood_Sj(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				write.table(x_complete, file = xConnection, row.names = FALSE, quote=F, append=TRUE)
				for(k in 1:K){
					myInd <- which(z == k)
					lMyInd <- length(myInd)
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection[[k]], append = TRUE)
					}
					cat(diag(SigmaINV.values[k,,]), " ", file = sigmainvConnection, append = TRUE)
					cat('\n', file = LambdaConnection[[k]], append = TRUE)
				}
				cat('\n', file = sigmainvConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(xConnection)
	for(k in 1:K){
#		close(regulonExpressionConnection[[k]])
		close(LambdaConnection[[k]])
	}
	setwd("../")

}




log_dirichlet_pdf <- function(alpha, weights){
	normConstant <- sum( lgamma(alpha) ) - lgamma( sum(alpha) )
	pdf <- sum( (alpha - 1)*log(weights) ) - normConstant
	return(pdf)
}

fabMix <- function(sameSigma = TRUE, dirPriorAlphas, rawData, outDir, Kmax, mCycles, burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize, thinning, zStart, nIterPerCycle, gibbs_z = 1){
	cat("          ____      __    __  ____     ", "\n")
	cat("         / __/___ _/ /_  /  |/  (_)  __", "\n")
	cat("        / /_/ __ `/ __ \\/ /|_/ / / |/_/", "\n")
	cat("       / __/ /_/ / /_/ / /  / / />  <  ", "\n")
	cat("      /_/  \\__,_/_.___/_/  /_/_/_/|_|  ", "\n")

	missingRowsIndex <- which(is.na(rowSums(rawData)) == TRUE)
	nMissingRows <- length( missingRowsIndex ) 
	if(nMissingRows > 0){
		stop("The data contains missing values. Use the `fabMix_missing_values()` function.")
	}
	dev.new(width=15, height=5)
	if(missing(Kmax)){Kmax <- 20}
	if(missing(nIterPerCycle)){nIterPerCycle = 10}
	if(missing(zStart)){zStart = FALSE}
	if( missing(dirPriorAlphas) ){
		nChains <- 8
		dN <- 1
		dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
#		dirPriorAlphas <- c((1/Kmax)*exp(seq(0,2,length = 10)),0.5,1)
#		dirPriorAlphas <- (1/Kmax)*exp( seq(0,4,length = 12) )
#		dirPriorAlphas <- seq(1,5,length = 8)/Kmax
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
	p <- dim(rawData)[2]
	n <- dim(rawData)[1]
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
	if(sameSigma == TRUE){
		cat(paste0('-    Parameterization: Same error variance per component'),'\n')
	}else{
		cat(paste0('-    Parameterization: Different error variance per component'),'\n')
	}
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
	d_per_cluster = 2*p + p*q + q*(q-1)/2
	initialAlphas <- seq(d_per_cluster/2, d_per_cluster, length = nChains)
	if(sameSigma == TRUE){
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 100, thinning = 1, burn = 99, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = 0.05)
		}
	}else{
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_Sj(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 100, thinning = 1, burn = 99, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = 0.05)
		}
	}
	cat(paste(' OK'),'\n')
	cat(paste('-    (2) Initializing the actual model from the previously obtained values... '))
	if(sameSigma == TRUE){
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 300, thinning = 1, burn = 299, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z)
		}
	}else{
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_Sj(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 300, thinning = 1, burn = 299, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z)
		}
	}
	cat(paste(' OK'),'\n')
	for(myChain in 1:nChains){
		kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
	}
	file_names <- list.files(outputDirs[1])	# the name of files are the same for each folder

	#	connections for saving the MCMC output corresponding to the target posterior distribution
	zConnection_target <- file(paste0(getwd(),"/zValues.txt"),open = "w")
	yConnection_target <- file(paste0(getwd(),"/yValues.txt"),open = "w")
	sigmainvConnection_target <- file(paste0(getwd(),"/sigmainvValues.txt"),open = "w")
	omegainvConnection_target <- file(paste0(getwd(),"/omegainvValues.txt"),open = "w")
	muConnection_target <- file(paste0(getwd(),"/muValues.txt"),open = "w")
	wConnection_target <- file(paste0(getwd(),"/wValues.txt"),open = "w")
	LambdaConnection_target <- vector('list',length = Kmax)
	cllConnection_target <- file(paste0(getwd(),"/cllValues.txt"),open = "w")
	for(k in 1:Kmax){
		LambdaConnection_target[[k]] <- file(paste0(getwd(),"/LambdaValues",k,".txt"),open = "w")  #K x p x q per iteration
	}




	cat(paste('-    (3) Running the sampler... '),'\n')
	#	main loops
	weights <- array(data = NA, dim = c(2, Kmax))
	#par(mfrow = c(1, 2))
	bb <- nIterPerCycle - 1
	for( iteration in 2:mCycles ){
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_Sj(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
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
			cat(paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.'), '\n')
		}
		if(iteration %% thinning == 0){
			if(iteration > burnCycles){
				y        <- as.numeric(read.table('alpha_1/yValues.txt'))
				w        <- as.numeric(read.table('alpha_1/wValues.txt')) 
				mu       <- as.numeric(read.table('alpha_1/muValues.txt'))
				omegainv <- as.numeric(read.table('alpha_1/omegainvValues.txt'))
				sigmainv <- as.numeric(read.table('alpha_1/sigmainvValues.txt')) 
				cll      <- as.numeric(read.table('alpha_1/k.and.logl.Values.txt') )[2]
				for(k in 1:Kmax){
					Lambda <- as.numeric( read.table(paste0('alpha_1/LambdaValues', k, '.txt') ) )
					cat(Lambda, file = LambdaConnection_target[[k]], '\n', append = TRUE)
				}
				cat(z       , file = zConnection_target, '\n', append = TRUE)
				cat(y       , file = yConnection_target, '\n', append = TRUE)
				cat(w       , file = wConnection_target, '\n', append = TRUE)
				cat(mu      , file = muConnection_target, '\n', append = TRUE)
				cat(omegainv, file = omegainvConnection_target, '\n', append = TRUE)
				cat(sigmainv, file = sigmainvConnection_target, '\n', append = TRUE)
				cat(cll     , file = cllConnection_target, '\n', append = TRUE)
			}
		}
	}
	stopImplicitCluster()
	close(zConnection_target)
	close(yConnection_target)
	close(wConnection_target)
	close(muConnection_target)
	close(sigmainvConnection_target)
	close(omegainvConnection_target)
	close(cllConnection_target)
	for(k in 1:Kmax){
		close(LambdaConnection_target[[k]])
	}
	keepedSeq <- seq(burnCycles + thinning, mCycles, by = thinning)
	if( burnCycles > 0){
		burnedSeq <- 1:burnCycles
		write.table(file = 'burn_in_period_k.txt', kValues[burnedSeq, ])
	}
	write.table(file = 'kValues.txt', kValues[keepedSeq, ], quote = FALSE, row.names = FALSE, col.names = dirPriorAlphas)
	setwd("../")
	cat('-    DONE.','\n')
}




fabMix_missing_values <- function(sameSigma = TRUE, dirPriorAlphas, rawData, outDir, Kmax, mCycles, burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize, thinning, zStart, nIterPerCycle, gibbs_z = 1){
	cat("          ____      __    __  ____     ", "\n")
	cat("         / __/___ _/ /_  /  |/  (_)  __", "\n")
	cat("        / /_/ __ `/ __ \\/ /|_/ / / |/_/", "\n")
	cat("       / __/ /_/ / /_/ / /  / / />  <  ", "\n")
	cat("      /_/  \\__,_/_.___/_/  /_/_/_/|_|  ", "\n")
	dev.new(width=15, height=5)
	if(missing(Kmax)){Kmax <- 20}
	if(missing(nIterPerCycle)){nIterPerCycle = 10}
	if(missing(zStart)){zStart = FALSE}
	if( missing(dirPriorAlphas) ){
		nChains <- 8
		dN <- 1
		dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
	}
	nChains <- length(dirPriorAlphas)
	if( range(diff(order(dirPriorAlphas)))[1] != 1){stop('dirPriorAlphas should be in increasing order.')}
	if( range(diff(order(dirPriorAlphas)))[2] != 1){stop('dirPriorAlphas should be in increasing order.')}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 1}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	p <- dim(rawData)[2]
	n <- dim(rawData)[1]
	dir.create(outDir)
	setwd(outDir)
	registerDoParallel(cores = nChains)
	outputDirs <- paste0('alpha_',1:nChains)
	originalX <- rawData
	x_data <- originalX
#	missing
	missingRowsIndex <- which(is.na(rowSums(x_data)) == TRUE)
	nMissingRows <- length( missingRowsIndex ) 
	missing_entries <- apply(x_data[missingRowsIndex, ], 1, function(y){which( is.na(y) == TRUE)})
	for(i in 1:nMissingRows){
		missing_entries[[i]] <- c(missingRowsIndex[i], missing_entries[[i]])
	}
	cat(paste0('-    Found ',nMissingRows,' rows containing missing values'),'\n')
#
	if( missing(thinning) ){thinning = 1}
	if( thinning < 1 ){ stop('thinning should be larger than or equal to 1.') }
	thinning <- floor(thinning)
	if( missing(normalize) ){normalize <- TRUE}
	cat(paste0("-    p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	if(sameSigma == TRUE){
		cat(paste0('-    Parameterization: Same error variance per component'),'\n')
	}else{
		cat(paste0('-    Parameterization: Different error variance per component'),'\n')
	}
	cat(paste0('-    Using Nchains = ', nChains),'\n')
	cat(paste0('-    Target posterior distribution corresponds to alpha = ', dirPriorAlphas[1]),'\n')
	if( normalize == TRUE ){
		x_data <- scale(originalX, center = T, scale = apply(originalX, 2, sd, na.rm = TRUE))
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
	d_per_cluster = 2*p + p*q + q*(q-1)/2
	initialAlphas <- seq(d_per_cluster/2, d_per_cluster, length = nChains)
	if(sameSigma == TRUE){
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 100, thinning = 1, burn = 99, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = 0.05)
		}
	}else{
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_Sj_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 100, thinning = 1, burn = 99, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = 0.05)
		}
	}
	cat(paste(' OK'),'\n')
	cat(paste('-    (2) Initializing the actual model from the previously obtained values... '))
	if(sameSigma == TRUE){
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 300, thinning = 1, burn = 299, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z)
		}
	}else{
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_Sj_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 300, thinning = 1, burn = 299, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z)
		}
	}
	cat(paste(' OK'),'\n')
	for(myChain in 1:nChains){
		kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
	}
	file_names <- list.files(outputDirs[1])	# the name of files are the same for each folder

	#	connections for saving the MCMC output corresponding to the target posterior distribution
	zConnection_target <- file(paste0(getwd(),"/zValues.txt"),open = "w")
	yConnection_target <- file(paste0(getwd(),"/yValues.txt"),open = "w")
	sigmainvConnection_target <- file(paste0(getwd(),"/sigmainvValues.txt"),open = "w")
	omegainvConnection_target <- file(paste0(getwd(),"/omegainvValues.txt"),open = "w")
	muConnection_target <- file(paste0(getwd(),"/muValues.txt"),open = "w")
	wConnection_target <- file(paste0(getwd(),"/wValues.txt"),open = "w")
	LambdaConnection_target <- vector('list',length = Kmax)
	cllConnection_target <- file(paste0(getwd(),"/cllValues.txt"),open = "w")
	for(k in 1:Kmax){
		LambdaConnection_target[[k]] <- file(paste0(getwd(),"/LambdaValues",k,".txt"),open = "w")  #K x p x q per iteration
	}




	cat(paste('-    (3) Running the sampler... '),'\n')
	#	main loops
	weights <- array(data = NA, dim = c(2, Kmax))
	#par(mfrow = c(1, 2))
	bb <- nIterPerCycle - 1
	xMean <- array(data = 0, dim = c(n, p))
	for( iteration in 2:mCycles ){
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_Sj_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
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
			cat(paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.'), '\n')
		}
		if(iteration %% thinning == 0){
			if(iteration > burnCycles){
				y        <- as.numeric(read.table('alpha_1/yValues.txt'))
				w        <- as.numeric(read.table('alpha_1/wValues.txt')) 
				mu       <- as.numeric(read.table('alpha_1/muValues.txt'))
				omegainv <- as.numeric(read.table('alpha_1/omegainvValues.txt'))
				sigmainv <- as.numeric(read.table('alpha_1/sigmainvValues.txt')) 
				cll      <- as.numeric(read.table('alpha_1/k.and.logl.Values.txt') )[2]
				xCurrent <- as.matrix(read.table('alpha_1/x_complete.txt', header = TRUE) )
				xMean    <- xMean + xCurrent
				for(k in 1:Kmax){
					Lambda <- as.numeric( read.table(paste0('alpha_1/LambdaValues', k, '.txt') ) )
					cat(Lambda, file = LambdaConnection_target[[k]], '\n', append = TRUE)
				}
				cat(z       , file = zConnection_target, '\n', append = TRUE)
				cat(y       , file = yConnection_target, '\n', append = TRUE)
				cat(w       , file = wConnection_target, '\n', append = TRUE)
				cat(mu      , file = muConnection_target, '\n', append = TRUE)
				cat(omegainv, file = omegainvConnection_target, '\n', append = TRUE)
				cat(sigmainv, file = sigmainvConnection_target, '\n', append = TRUE)
				cat(cll     , file = cllConnection_target, '\n', append = TRUE)
			}
		}
	}
	stopImplicitCluster()
	close(zConnection_target)
	close(yConnection_target)
	close(wConnection_target)
	close(muConnection_target)
	close(sigmainvConnection_target)
	close(omegainvConnection_target)
	close(cllConnection_target)
	for(k in 1:Kmax){
		close(LambdaConnection_target[[k]])
	}
	keepedSeq <- seq(burnCycles + thinning, mCycles, by = thinning)
	if( burnCycles > 0){
		burnedSeq <- 1:burnCycles
		write.table(file = 'burn_in_period_k.txt', kValues[burnedSeq, ])
	}
	write.table(file = 'kValues.txt', kValues[keepedSeq, ], quote = FALSE, row.names = FALSE, col.names = dirPriorAlphas)
	write.table(file = 'xMeanEstimate.txt', xMean/(mCycles - burnCycles), quote = FALSE, row.names = FALSE)
	setwd("../")
	cat('-    DONE.','\n')
}




observed.log.likelihood0 <- function(x_data, w, mu, Lambda, Sigma, z){
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	ct <- -(p/2)*log(2*pi)

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
                diag(x_var) <- diag(x_var) + Sigma
#                x_var <- try(solve(x_var), TRUE)
#                loggedValues[ ,as.character(k)] <- log(newW[k]) - 0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) + ct
                loggedValues[ ,as.character(k)] <- log(newW[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
        }
        lMax <- apply(loggedValues, 1, max)
        if( length(alive) == 1 ){
                logL <- sum(lMax + log( exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ) ) )
        }else{
                logL <- sum(lMax + log( rowSums(exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ))))
        }
        return( logL )
}


observed.log.likelihood0_Sj <- function(x_data, w, mu, Lambda, Sigma, z){
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	ct <- -(p/2)*log(2*pi)

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
#                x_var <- try(solve(x_var), TRUE)
#                loggedValues[ ,as.character(k)] <- log(newW[k]) - 0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) + ct
                loggedValues[ ,as.character(k)] <- log(newW[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
        }
        lMax <- apply(loggedValues, 1, max)
        if( length(alive) == 1 ){
                logL <- sum(lMax + log( exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ) ) )
        }else{
                logL <- sum(lMax + log( rowSums(exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ))))
        }
        return( logL )
}





getStuffForDIC <- function(sameSigma = TRUE, x_data, outputFolder, q, burn, Km, normalize, discardLower){
	cat(paste0('-    (4) Computing information criteria for q = ', q), '\n')
	if(missing(normalize)){normalize = TRUE}
	if(normalize){
		x_data <- scale(x_data, center = TRUE, scale = TRUE)
		cat('-    NOTE: using standardizing data.','\n')
	}
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]

	if(missing(Km)){Km <- 20}
	if(missing(burn)){burn <- 0}
	setwd(outputFolder)
	cat(paste0('         - Entering directory: ', getwd()),'\n')
	z <- as.matrix(read.table("zValues.txt"))
	logl <- read.table("cllValues.txt")
        tmp  <- apply(z,1,function(y){length(table(y))})
        logl <- cbind(tmp, logl)
	if(burn > 0){
		z <- z[-(1:burn),]
		logl <- logl[-(1:burn),]
	}
	cat(paste0('            Nclusters:    ', paste(as.character(names(table(logl[,1]))), collapse="     ") ), '\n')
	cat(paste0('            Frequency:    ', paste(as.character(as.numeric(table(logl[,1]))), collapse="    ") ), '\n')
#	K <- Km
#	kSelected <- K
#	index <- 1:dim(z)[1]
#	Kindex <- index
        K <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
        kSelected <- K
        index <- which(logl[,1] == K)
        Kindex <- index
	m <- length(index)

	#this is artificial

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
	mcmc <- array(data = NA, dim = c(m,Km,J))
	for(k in 1:Km){
	#	lMean <- apply(l,2,mean)
	#	lMean.matrix <- matrix(lMean,nrow = p, ncol = q, byrow=TRUE)
		l <- as.matrix(read.table(paste0("LambdaValues",k,".txt")))
#		if(q == 1){ l <- array(l, dim = c(length(l) , 1))}
		if(burn > 0){
		l <- l[-(1:burn),]}
		mcmc[,k,] <- l[index,]
	}
	lambda.perm.mcmc <- permute.mcmc(mcmc, ls$permutations$ECR)$output
	for(k in 1:Km){
		lMean <- apply(lambda.perm.mcmc[,k,],2,mean)
	}
	mu <- read.table("muValues.txt")# auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
	if(burn > 0){
		mu <- mu[-(1:burn),] 
	}
	mu <- mu[Kindex,]
	mu.mcmc <- array(data = NA, dim = c(m,Km,p))
	for(k in 1:Km){
		mu.mcmc[,k,] <- as.matrix(mu[,k + Km*((1:p)-1)])
	}
	mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations$ECR)$output
	#
	if(sameSigma == TRUE){
		SigmaINV <- as.matrix(read.table("sigmainvValues.txt"))# auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
		if(burn > 0){
			SigmaINV <- SigmaINV[-(1:burn),] 
		}
		SigmaINV <- SigmaINV[Kindex, ] 
		SigmaINV.mcmc <- SigmaINV
	}else{
		SigmaINV <- as.matrix(read.table("sigmainvValues.txt")) # auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
		if(burn > 0){
			SigmaINV <- SigmaINV[-(1:burn),] 
		}
		SigmaINV <- SigmaINV[Kindex,  ] 
		SigmaINV.mcmc <- array(data = NA, dim = c(m,Km,p))
		for(k in 1:Km){
		        SigmaINV.mcmc[,k,] <- as.matrix(SigmaINV[,((k-1)*p + 1):(k*p)])
		}
		SigmaINV.mcmc <- permute.mcmc(SigmaINV.mcmc, ls$permutations$ECR)$output

	}
	Sigma.mcmc <- 1/SigmaINV.mcmc

	#SigmaINV.mean <- as.numeric(apply(SigmaINV,2,mean))
	w.mcmc <- as.matrix(read.table("wValues.txt"))
	w.mcmc <- array(w.mcmc, dim = c(dim(w.mcmc)[1], Km, 1))
	if(burn > 0){
		w.mcmc <- w.mcmc[-(1:burn),,]
		w.mcmc <- w.mcmc[Kindex,]
	}else{
		w.mcmc <- w.mcmc[Kindex,,]
	}
	w.mcmc <- array(w.mcmc[,1:Km],dim = c(length(Kindex),Km,1))
	w.mcmc <- permute.mcmc(w.mcmc, ls$permutations$"ECR")$output

	lValues <- numeric(m)
	cll <- 0
	aic <- 0
	bic <- 0
	maxL <- 1
	i <- 1
	lambda.current <- array(data = NA, dim = c(Km,p,q))
	if(sameSigma == TRUE){
		Sigma.current <- Sigma.mcmc[i, ]
	}else{
	        Sigma.current <- Sigma.mcmc[i, , ]		
	}
	mu.current <- mu.mcmc[i,,]
	if( Km == 1 ){ 	
		Sigma.current <- array(Sigma.current, dim = c(1, p))
		mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
	}
	for(k in 1:Km){
		ldraw <- lambda.perm.mcmc[i,k,]
		lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
		for(i1 in 1:p){
			if( Sigma.current[ i1] > 100 ){ cat(paste0('oops: ', i),'\n'); Sigma.current[i1] <- 100 }
		}
	}
	if(sameSigma == TRUE){
		obsL <- observed.log.likelihood0(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
	}else{
		obsL <- observed.log.likelihood0_Sj(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
	}
	lValues[i] <- obsL
	maxL <- obsL
	cll <- cll + obsL
	aic <- aic + obsL
	bic <- bic + obsL
	iterMax <- i
	if(sameSigma == TRUE){
		for(i in 2:m){
	#		cat(paste0("i  = ", i), "\n")
			lambda.current <- array(data = NA, dim = c(Km,p,q))
			Sigma.current <- Sigma.mcmc[i, ]
			mu.current <- mu.mcmc[i,,]
			if( Km == 1 ){ 	
				Sigma.current <- array(Sigma.current, dim = c(1, p))
				mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
			}
			for(k in 1:Km){
	#			cat(paste0("  k  = ", k), "\n")
				ldraw <- lambda.perm.mcmc[i,k,]
				lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
				for(i1 in 1:p){
					if( Sigma.current[ i1] > 100 ){ cat(paste0('oops: ', i),'\n'); Sigma.current[ i1] <- 100 }
				}
			}
			obsL <- observed.log.likelihood0(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
			lValues[i] <- obsL
			if( obsL > maxL ){
				maxL <- obsL	
				iterMax <- i
			}
			cll <- cll + obsL
			aic <- aic + obsL
			bic <- bic + obsL
		}
	}else{
		for(i in 2:m){
	#		cat(paste0("i  = ", i), "\n")
			lambda.current <- array(data = NA, dim = c(Km,p,q))
			Sigma.current <- Sigma.mcmc[i, , ]
			mu.current <- mu.mcmc[i,,]
			if( Km == 1 ){ 	
				Sigma.current <- array(Sigma.current, dim = c(1, p))
				mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
			}
			for(k in 1:Km){
	#			cat(paste0("  k  = ", k), "\n")
				ldraw <- lambda.perm.mcmc[i,k,]
				lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
				for(i1 in 1:p){
					if( Sigma.current[k, i1] > 100 ){ Sigma.current[k, i1] <- 100 }
				}
			}
			obsL <- observed.log.likelihood0_Sj(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
			lValues[i] <- obsL
			if( obsL > maxL ){
				maxL <- obsL	
				iterMax <- i
			}
			cll <- cll + obsL
			aic <- aic + obsL
			bic <- bic + obsL
		}
	}
	if(missing(discardLower)){ discardLower <- 0.1 }
	if ( discardLower == FALSE){
		cll <- cll/m
	}else{
		cll <- mean( lValues[which( lValues > as.numeric(quantile(lValues, discardLower)) )] )
	}
#	dic_classic <- -4*cll + 2*logL.theta_hat
#	dic_star <- -6*cll + 4*logL.theta_hat
#	dic_classicMAP <- -4*cll + 2*logL.theta_map
#	dic_starMAP <- -6*cll + 4*logL.theta_map
	dic_classicMAP <- -4*cll + 2*maxL
	dic_starMAP <- -6*cll + 4*maxL
	dic_classic <- dic_classicMAP
	dic_star <- dic_starMAP
	if(sameSigma == TRUE){
		aic <- -2*aic/m + 2*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 ) 
		bic <- -2*bic/m + log(n)*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 )
		aic_MAX <- -2*maxL + 2*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 ) 
		bic_MAX <- -2*maxL + log(n)*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 )
	}else{
		aic <- -2*aic/m + 2*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 ) 
		bic <- -2*bic/m + log(n)*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 )
		aic_MAX <- -2*maxL + 2*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 ) 
		bic_MAX <- -2*maxL + log(n)*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 )
	}
	dic <- c(aic, bic, dic_classic, dic_star, dic_classicMAP, dic_starMAP, aic_MAX, bic_MAX)
	names(dic) <- c('AIC', 'BIC', 'DIC1', 'DIC*2', 'DIC', 'DIC_2', 'AIC_map', 'BIC_map')	
#	write.table(file = 'informationCriteria_map_model.txt', dic[c(1,2,5,6,7,8)], col.names = paste0('q_',q), quote = FALSE)
	write.table(file = 'informationCriteria_map_model.txt', dic[c(5,6,7,8)], col.names = paste0('q_',q), quote = FALSE)
	write.table(file = 'lValues_map.txt', lValues, quote = FALSE)
	setwd("../")
	cat(paste0('         - Information criteria written to `', outputFolder,'/informationCriteria_map_model.txt`.'), '\n')
}



dealWithLabelSwitching <- function(sameSigma = TRUE, x_data, outputFolder, q, burn, z.true, compute_regularized_expression, Km){
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	if(missing(Km)){Km <- 20}
	if(missing(burn)){burn <- 0}
	if(missing(compute_regularized_expression)){ compute_regularized_expression = FALSE }
	cat(paste0('-    (5) Dealing with label switching for q = ', q), '\n')
	setwd(outputFolder)
	cat(paste0('         * Entering directory: ', getwd()),'\n')
	z <- as.matrix(read.table("zValues.txt"))
	logl <- read.table("kValues.txt", header=T)
	if(burn > 0){
		logl <- logl[-(1:burn), ]
		z <- z[-(1:burn),]
	}
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
		skiniko <- FALSE
		if(K < Km){
			skiniko <- TRUE
			kLAST <- K 
			newZ <- t(apply(z, 1, function(y){ myI <- which(y == K); y[myI] <- rep(Km, length(myI));return(y) } ) )
			z <- newZ
			newZ <- zPivot
			myI <- which(newZ == K)
			zPivot[myI] <- rep(Km, length(myI))
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
			if(burn > 0){
				l <- l[-(1:burn),]
			}
			mcmc[,k,] <- l[index,]
		}
		if(skiniko == TRUE){
			tmp1 <- mcmc[,kLAST,]
			tmp2 <-  mcmc[,Km,]	
			mcmc[,kLAST,] <- tmp2
			mcmc[,Km,] <- tmp1
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
		mu <- read.table("muValues.txt")# auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
		if(burn > 0){
			mu <- mu[-(1:burn),] 
		}
		mu <- mu[Kindex,]

		mu.mcmc <- array(data = NA, dim = c(m,K,p))
		for(k in 1:K){
			mu.mcmc[,k,] <- as.matrix(mu[,k + K*((1:p)-1)])
		}
		if(skiniko == TRUE){
			tmp1 <- mu.mcmc[,kLAST,]
			tmp2 <-  mu.mcmc[,Km,]	
			mu.mcmc[,kLAST,] <- tmp2
			mu.mcmc[,Km,] <- tmp1
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
		if(sameSigma == TRUE){
			sigmaINV <- as.matrix(read.table("sigmainvValues.txt"))# auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
			if(burn > 0){
				sigmaINV <- sigmaINV[-(1:burn),] 
			}
			sigmaINV <- sigmaINV[Kindex, ] 
		}else{	
		
		        SigmaINV <- as.matrix(read.table("sigmainvValues.txt"))
			if(burn > 0){
				SigmaINV <- SigmaINV[-(1:burn),] 
			}
			SigmaINV <- SigmaINV[Kindex,]

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
		}
#		w.mcmc <- array(as.matrix(read.table("wValues.txt")[max(burn) + Kindex, ]),dim = c(length(Kindex),K,1))

		w.mcmc <- as.matrix(read.table("wValues.txt"))
		w.mcmc <- array(w.mcmc, dim = c(dim(w.mcmc)[1], K, 1))
		if(burn > 0){
			w.mcmc <- w.mcmc[-(1:burn),,]
			w.mcmc <- w.mcmc[Kindex,]
		}else{
			w.mcmc <- w.mcmc[Kindex,,]
		}
		w.mcmc <- array(w.mcmc[,1:K],dim = c(length(Kindex),K,1))
		w.mcmc_raw <- w.mcmc
		if(skiniko == TRUE){
			tmp1 <- w.mcmc[,kLAST,]
			tmp2 <-  w.mcmc[,Km,]	
			w.mcmc[,kLAST,] <- tmp2
			w.mcmc[,Km,] <- tmp1
		}
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
		zMAP <- as.matrix(read.table("singleBestClusterings.txt",header=TRUE))[1, ]


		mapAllocationsPosteriorProbs <- numeric(n)
		for(i in 1:n){
			mapAllocationsPosteriorProbs[i] <- length(which(allocationsECR[,i] == zMAP[i]))
		}
		mapAllocationsPosteriorProbs <- mapAllocationsPosteriorProbs/dim(allocationsECR)[1]
		write.table(mapAllocationsPosteriorProbs, file = "classificationProbabilities.txt")
		cat(paste0('         * write file: `classificationProbabilities.txt`'),'\n')

		aliveClusters <- as.numeric(names(table(ls$clusters[1,])))

		covmat <- array(data = 0, dim = c( length(aliveClusters), p, p ))
		rownames(covmat) <- as.character(aliveClusters)
		if(sameSigma == TRUE){
			for(iter in 1:m){
				for(k in aliveClusters){
					lambda_k <- matrix( lambda.perm.mcmc[iter, k, ], nrow = p, ncol = q, byrow=T)
					covmat[as.character(k), , ] <- covmat[as.character(k), , ] + lambda_k %*% t(lambda_k) 
					diag(covmat[as.character(k), , ]) <- diag(covmat[as.character(k), , ]) + as.numeric(1/sigmaINV[iter, ])
				}
			}
		}else{
			for(iter in 1:m){
				for(k in aliveClusters){
					lambda_k <- matrix( lambda.perm.mcmc[iter, k, ], nrow = p, ncol = q, byrow=T)
					covmat[as.character(k), , ] <- covmat[as.character(k), , ] + lambda_k %*% t(lambda_k) 
					diag(covmat[as.character(k), , ]) <- diag(covmat[as.character(k), , ]) + as.numeric(Sigma.mcmc[iter,k,])
				}
			}
		}
		covmat <- covmat/m
		for(k in 1:length(aliveClusters)){
			write.table(covmat[k, , ], file = paste0("estimated_cov_cluster_",k,".txt"))
			cat(paste0('         * write file: `estimated_cov_cluster_',k,'.txt`'),'\n')
			write.table(lambda.map[k, , ], file = paste0('lambda_map_',k,'.txt'))
			cat(paste0('         * write file: `lambda_map_',k,'.txt`'),'\n')
		}




		if( compute_regularized_expression == TRUE ){
			cat(paste("-    computing regularized expressions..."),'\n')
			yValues <- read.table("yValues.txt")
			if(burn > 0){
				yValues <- yValues[-(1:burn), ]
			}
			yValues <- yValues[Kindex, ]
			regularizedExpression <- array(data = 0, dim = c(length(aliveClusters), p))
			regularizedExpression2 <- array(data = 0, dim = c(length(aliveClusters), p, q))			
			rownames(regularizedExpression) <- rownames(regularizedExpression2) <- as.character(aliveClusters)
			for(iter in 1:m){
				for(k in aliveClusters){
					lambda_k <- matrix( lambda.perm.mcmc[iter, k, ], nrow = p, ncol = q, byrow=T)
					yMean <- array(data = 0, dim = c(q,1))
					index_k <- as.numeric(which( allocationsECR[iter, ] == k ))
					for(i in index_k ){
						yMean <- yMean + t(yValues[iter , (1:q - 1)*n + i])
					}
					yMean <- yMean/length(index_k)
					tmp <- t(apply(lambda_k,1, function(y){ y*yMean }))
					regularizedExpression2[as.character(k), , ] <- regularizedExpression2[as.character(k), , ] + tmp
					regularizedExpression[as.character(k), ] <- regularizedExpression[as.character(k), ] + rowSums(tmp)
				}
				if(iter %% 50 == 0){cat(paste('iter = ', iter),'\n')}
			}
			regularizedExpression <- regularizedExpression/m
			regularizedExpression2 <- regularizedExpression2/m
			for(j in 1:q){
				write.table(
				regularizedExpression2[,,j], 
				file = paste0("estimated_regularized_expression_per_cluster_",j,".txt"))
				cat(paste0('         * write file: `estimated_regularized_expression_per_cluster_',j,'.txt`'),'\n')

			}
			write.table(regularizedExpression, file = "estimated_regularized_expression_per_cluster.txt")
			cat(paste0('         * write file: `estimated_regularized_expression_per_cluster.txt`'),'\n')
		}


		cat(paste0('-    Done.'), '\n')
	}
	setwd("../")

}



simData <- function(sameSigma = TRUE, p, q, K.true, n, loading_means, loading_sd, sINV_values){
	if(missing(p)){p = 40}
	if(missing(q)){q = 4}
	if(missing(K.true)){K.true = 6}
	if(missing(n)){n = 200}
	if(missing(sINV_values)){
		if(sameSigma){ 
			sINV_values = rgamma(p, shape = 1, rate = 1) 
		}else{
			sINV_values = matrix(rgamma(K.true*p, shape = 1, rate = 1), nrow = K.true, ncol = p )
		}
	}
	if( missing(loading_means) ){ loading_means <- c(-30,-20,-10,10, 20, 30) }
	if( missing(loading_sd) ){ loading_sd <- rep(2, length(loading_means)) }
	if ( length(loading_means) !=  length(loading_sd) ){
		stop("`loading_means` and `loading_sd` should have same length.")
	}
	cat(paste0("Simulation parameters:"),'\n')
	if(q >= p){stop("q should not be greater than p")}
	cat(paste0("   n = ", n),'\n')
	cat(paste0("   p = ", p),'\n')
	cat(paste0("   q = ", q),'\n')
	cat(paste0("   K = ", K.true),'\n')
	w.true <- myDirichlet(rep(10,K.true))
	z.true <- sample(K.true,n,replace = TRUE, prob = w.true)
	Lambda.true <- array(data = rnorm(K.true*p*q, mean = 0, sd = 1), dim = c(K.true,p,q))
	mu.true <- array(data = 0, dim = c(K.true,p))
	for(k in 1:K.true){
		u <- runif(1)
		subROW <- floor(p/q)
		for(j in 1:q){
			meanCOL <- rep(0,p)
			pickAnIndex <- sample(length(loading_means), 1)
			meanCOL[ (j-1)*subROW + 1:subROW] <- loading_means[pickAnIndex]
			Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = loading_sd[pickAnIndex] )
			if(j > 1)  {
				Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)
			}
		}
		u <- runif(1)
		if(u < 1/3){
			mu.true[k, ] <- 20*sin(seq(0,k*pi, length = p))
		}else{
			if(u < 2/3){
				mu.true[k, ] <- 20*cos(seq(0,k*pi, length = p))
			}else{
				mu.true[k, ] <- 40*(sin(seq(0,k*pi, length = p)))^2 - 40*(cos(seq(0,k*pi, length = p)))^2
			}
		}
	}


	if(sameSigma == TRUE){
		SigmaINV.true <- array(data = 0, dim = c(p,p))
		diag(SigmaINV.true) <- sINV_values
		Sigma.true <- SigmaINV.true
		diag(Sigma.true) <- 1/diag(SigmaINV.true)

	}else{
		SigmaINV.true <- array(data = 0, dim = c(K.true, p,p))
		Sigma.true <- SigmaINV.true
		for(k in 1:K.true){
			diag(SigmaINV.true[k,,]) <- sINV_values[k,]
			diag(Sigma.true[k,,]) <- 1/diag(SigmaINV.true[k,,])
		}
	}
	y.true <- array(data = 0, dim = c(n,q))
	x_data <- array(data = 0, dim = c(n,p))
	ly <- q
	for(i in 1:n){
		y.true[i,] <- rnorm(ly,mean = 0,sd = 1)
		j <- z.true[i]
		if(q == 1){
			x_mean <- mu.true[j,] + Lambda.true[j, , ] %*% array(y.true[i, ], dim = c(q,1))
		}else{
			x_mean <- mu.true[j,] + Lambda.true[j, , ] %*% y.true[i, ]
		}
		if(sameSigma){
			x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true)
		}else{
			x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true[j,,])
		}
	}
	matplot(t(x_data), type = "l", col = z.true, lty = 1)
	legend("bottomleft", paste0("cluster ",1:K.true, ": ",as.character(as.numeric(table(z.true)))), col = 1:K.true, lty = 1)
	results <- vector('list', length = 7)
	results[[1]] <- x_data
	results[[2]] <- z.true
	results[[3]] <- Lambda.true
	results[[4]] <- mu.true 
	results[[5]] <- Sigma.true
	results[[6]] <- y.true
	results[[7]] <- w.true
	names(results) <- c("data", "class", "factorLoadings", "means", "variance","factors","weights")
	return(results)
}


