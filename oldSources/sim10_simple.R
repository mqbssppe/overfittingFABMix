library(MASS)
library("ellipse")
library('mclust')
library('flexmix')

# Simulate data from the model



# dirichlet function
myDirichlet <- function(alpha){
        k <- length(alpha)
        theta <- rgamma(k, shape = alpha, rate = 1)
        return(theta/sum(theta))
}

#set.seed(100)
K.true <- 10 
outFile <- file(paste0("nClusters_",K.true,"_99.txt"), "w")  # open an output file connection
 
#for(dIter in 1:10){
		dIter <- 99
		set.seed(454)
  		w.true <- myDirichlet(rep(100,K.true))
		w.true <- w.true/sum(w.true)
		p <- 40
		q <- 4
		n <- 1000
		z.true <- sample(K.true,n,replace = TRUE, prob = w.true)
		v_r <- numeric(p) #indicates the non-zero values of Lambdas
		Lambda.true <- array(data = 0, dim = c(K.true,p,q))
		mu.true <- array(data = 0, dim = c(K.true,p))
		mu.0 <- 30*sin((1:p)*pi/10)
		for(k in 1:floor(K.true/2 + 0.5)){
			# auto einai kai kala gia na min einai 0
			for( r in 1:p ){
				v_r[r] <- min(r,q)
				Lambda.true[k,r,1:v_r[r]] <- rnorm(v_r[r], mean = 0, sd = 0.1)  # -k + k * runif(v_r[r])
			}

			j <- 1
			meanCOL <- rep(0,p)
			meanCOL[1:10] <- 10
			Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  

			j <- 2
			meanCOL <- rep(0,p)
			meanCOL[11:20] <- 20
			Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  
			Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)

			j <- 3
			meanCOL <- rep(0,p)
			meanCOL[21:30] <- 10
			Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  
			Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)

			j <- 4
			meanCOL <- rep(0,p)
			meanCOL[31:40] <- 20
			Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  
			Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)
			mu.true[k, ] <- 100*k

		}


		for(k in (1+floor(K.true/2 + 0.5)):K.true){
			# auto einai kai kala gia na min einai 0
			for( r in 1:p ){
				v_r[r] <- min(r,q)
				Lambda.true[k,r,1:v_r[r]] <- rnorm(v_r[r], mean = 0, sd = 0.1)  # -k + k * runif(v_r[r])
			}

			j <- 1
			meanCOL <- rep(0,p)
			meanCOL[1:10] <- -40
			Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  

			j <- 2
			meanCOL <- rep(0,p)
			meanCOL[11:20] <- 10
			Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  
			Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)

			j <- 3
			meanCOL <- rep(0,p)
			meanCOL[21:30] <- 30
			Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  
			Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)

			j <- 4
			meanCOL <- rep(0,p)
			meanCOL[31:40] <- -20
			Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  
			Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)
			mu.true[k, ] <- 100*k

		}








		SigmaINV.true <- array(data = 0, dim = c(p,p))
		diag(SigmaINV.true) <- rgamma(p, shape = 5, rate = 5)
		#diag(SigmaINV.true) <- 1/((1:p))
		y.true <- array(data = 0, dim = c(n,q))
		Sigma.true <- SigmaINV.true
		diag(Sigma.true) <- 1/diag(SigmaINV.true)
		x_data <- array(data = 0, dim = c(n,p))
		for(i in 1:n){
			ly <- q
			y.true[i,] <- rnorm(ly,mean = 0,sd = 1)
			#y.true[i,] <- rnorm(ly,mean = 0,sd = sample(c(1,2,3,4),1))
			j <- z.true[i]
			x_mean <- mu.true[j,] + Lambda.true[j, , ] %*% y.true[i, ]
			x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true)
		}

		#ss <- array(data = 0, dim = c(p,p))
		#for(r in 1:p){
		#       ss[r,] <- rnorm(p)
		#}
		#ss <- ss%*%t(ss)
		#x_data <- mvrnorm(n = n, mu = rep(0,p), Sigma = ss)
		originalX <- x_data
		x_data <- scale(x_data, center = TRUE, scale = TRUE)
		par(mfrow = c(1,2))
		matplot( t(x_data),col = z.true, type = "l" )
		matplot( t(originalX),col = z.true, type = "l" )
	




