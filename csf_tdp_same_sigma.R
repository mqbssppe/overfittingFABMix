# computation of log{ f(d_il = 0| ...)/f(d_il = 1| ...) }
library('MASS')
library('mvtnorm')
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
compute_A_B_G_D_and_simulate_mu_Lambda <- function(SigmaINV, suff_statistics, OmegaINV, K, priorConst1, T_INV){
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

### compute only the matrices

compute_A_B_G_D <- function(SigmaINV, suff_statistics, OmegaINV,K, T_INV){
	A <- array(data = 0, dim = c(K,p,p))
	B <- mu <- array(data = 0, dim = c(K,p))
	G <- array(data = 0, dim = c(K,p,q,q))
	D <- array(data = 0, dim = c(K,p,q))
	#cat(paste(v_r), '\n')
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
		}
		#B[k,] <- SigmaINV %*% (suff_statistics$sx[k,] - Lambdas[k,,] %*% suff_statistics$sy[k,] ) + priorConst1
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
		center_x <- x_data - t(apply(y,1,function(tmp){return(mu[k,] + array(Lambda[k,,], dim = c(p,q)) %*% tmp)}))
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




#simulate z and mixture weights from standard Gibbs
update_z_b <- function(w, mu, Lambda, y, SigmaINV, K){
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



# using dmvnorm from package mvtnorm
update_z4 <- function(w, mu, Lambda, SigmaINV, K){
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



# collapsed over y
update_z2 <- function(w, mu, Lambda, SigmaINV, K){
	probs <- array(data = 0, dim =c(n,K))
	for(k in 1:K){
		center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
		x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
		diag(x_var) <- diag(x_var) + 1/diag(SigmaINV)
		#x_var <- solve(x_var)
		x_var <- try(solve(x_var), TRUE)
		if(is.numeric(x_var) == TRUE){
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


update_z3 <- function(w, mu, Lambda, y, SigmaINV, K){
	probs <- array(data = 0, dim =c(n,K))
	index <- 1:n
	rn <- floor((n-2)*runif(1)) + 2 #oxi n gia na min midenizetai to allo
	index2 <- sort(sample(n,rn,replace=FALSE))
	index1 <- setdiff(1:n,index2)
	for(k in 1:K){
		center_x <- x_data[index1,] - matrix(mu[k,], nrow = length(index1), ncol = p, byrow=TRUE)
		x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
		diag(x_var) <- diag(x_var) + 1/diag(SigmaINV)
		#x_var <- solve(x_var)
		x_var <- try(solve(x_var), TRUE)
		if(is.numeric(x_var) == TRUE){
			probs[index1,k] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
		}
		center_x <- x_data[index2,] - t(apply(array(y[index2,],dim = c(length(index2),q)),1,function(tmp){return(mu[k,] + array(Lambda[k,,],dim=c(p,q)) %*% tmp)}))
		probs[index2,k] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% SigmaINV %*% tmp) )}) + 0.5*sum(log(diag(SigmaINV))) 

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



update_SigmaINV <- function(z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
        SigmaINV <- array(data = 0, dim = c(p,p))
        s <- numeric(p)
        for(i in 1:n){
                tmp <- x_data[i, ] - (mu[z[i],] + array(Lambda[ z[i], , ], dim = c(p,q)) %*% y[i,])
                # apo ton pinaka tmp x tmp xreiazomai mono ta diagonia stoixeia
                #s <- s + tmp %*% t(tmp)
                s <- s + tmp^2
        }
        diag(SigmaINV) <- rgamma(p,shape = alpha_sigma + n/2, rate = beta_sigma + s/2)
        return(SigmaINV)
}


update_SigmaINV_faster <- function(z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
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

truncatedDirichletProcess <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, progressGraphs, start_values, q, zStart, gibbs_z){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
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
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
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
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z4(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K)
		}else{
			f2 <- update_z_b(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = SigmaINV.values, K = K)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y(x = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster(z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		#print(OmegaINV.constant)
#		f2 <- update_z4(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K)
#		f2 <- update_z_b(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
#					SigmaINV = SigmaINV.values, K = K)

#		u_v <- runif(1)
#		if(u_v < 0.5){
#			f2 <- update_z_b(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
#					SigmaINV = SigmaINV.values, K = K)
#		}else{
#			if(u < 1.95){

#		}
#			}else{
#			f2 <- update_z3(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
#						SigmaINV =array(SigmaINV.values,dim = c(K,p,p)), K = K)}
#		}

		if(iter %% thinning == 0){
			#cat(paste0(iter," iterations done."),"\n")
			#print(table(z))
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood(w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
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

library('doParallel')


log_dirichlet_pdf <- function(alpha, weights){
	normConstant <- sum( lgamma(alpha) ) - lgamma( sum(alpha) )
	pdf <- sum( (alpha - 1)*log(weights) ) - normConstant
	return(pdf)
}

heated_chains <- function(dirPriorAlphas, rawData, outDir, Kmax, mCycles, burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize, thinning, zStart, nIterPerCycle){
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
	d_per_cluster = 2*p + p*q + q*(q-1)/2
	initialAlphas <- seq(d_per_cluster/2, d_per_cluster, length = nChains)
	foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		truncatedDirichletProcess(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
			Kmax = Kmax, m = 100, thinning = 1, burn = 99, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
			alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, progressGraphs = FALSE, start_values = FALSE, gibbs_z = 0.05)
	}
	cat(paste(' OK'),'\n')
	cat(paste('-    (2) Initializing the actual model from the previously obtained values... '))
	foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		truncatedDirichletProcess(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
			Kmax = Kmax, m = 300, thinning = 1, burn = 299, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
			alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, progressGraphs = FALSE, start_values = TRUE, gibbs_z = 0.05)
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
		
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			truncatedDirichletProcess(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
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
			cat(paste0('-        iteration: ',iteration,'. Metropolis-Hastings acceptance rate: ', ar), '\n')
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



