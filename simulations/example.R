library(fabMix)
library(MASS)

dIter <- 1
set.seed(101)
q <- q.true <- 4  
K.true <- 6

w.true <- myDirichlet(rep(10,K.true))
w.true <- w.true/sum(w.true)
p <- 40
q <- 4
n <- 200
z.true <- sample(K.true,n,replace = TRUE, prob = w.true)
v_r <- numeric(p) #indicates the non-zero values of Lambdas
Lambda.true <- array(data = 0, dim = c(K.true,p,q))
mu.true <- array(data = 0, dim = c(K.true,p))
for(k in 1:K.true){
	# auto einai kai kala gia na min einai 0
	for( r in 1:p ){
		v_r[r] <- min(r,q)
		Lambda.true[k,r,1:v_r[r]] <- rnorm(v_r[r], mean = 0, sd = 1)  # -k + k * runif(v_r[r])
	}
	u <- runif(1)
	if(u < 0.5){
		j <- 1
		meanCOL <- rep(0,p)
		meanCOL[1:10] <- 10
		Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  

		j <- 2
		meanCOL <- rep(0,p)
		meanCOL[11:20] <- 10
		Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  
		Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)

		j <- 3
		meanCOL <- rep(0,p)
		meanCOL[21:30] <- 10
		Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  
		Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)

		j <- 4
		meanCOL <- rep(0,p)
		meanCOL[31:40] <- 10
		Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 1)  
		Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)
	}else{

		j <- 1
		meanCOL <- rep(0,p)
		meanCOL[1:10] <- -10
		Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 2)  

		j <- 2
		meanCOL <- rep(0,p)
		meanCOL[11:20] <- 20
		Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 2)  
		Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)

		j <- 3
		meanCOL <- rep(0,p)
		meanCOL[21:30] <- 20
		Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 2)  
		Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)

		j <- 4
		meanCOL <- rep(0,p)
		meanCOL[31:40] <- -10
		Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = 2)  
		Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)


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



SigmaINV.true <- array(data = 0, dim = c(p,p))
diag(SigmaINV.true) <- rgamma(p, shape = 0.1, rate = 0.1)
diag(SigmaINV.true) <- 1/((1:p))
y.true <- array(data = 0, dim = c(n,q))
Sigma.true <- SigmaINV.true
diag(Sigma.true) <- 1/diag(SigmaINV.true)
x_data <- array(data = 0, dim = c(n,p))
for(i in 1:n){
	ly <- q
	y.true[i,] <- rnorm(ly,mean = 0,sd = 1)
	j <- z.true[i]
	x_mean <- mu.true[j,] + Lambda.true[j, , ] %*% y.true[i, ]
	x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true)
}
originalX <- x_data

write.table(originalX, file = paste0('data_',K.true,'_',dIter,'.txt'))
write.table(z.true, file = paste0('z_',K.true,'_',dIter,'.txt'))
matplot(t(originalX), type = "l", col = z.true)
table(z.true)


#}







simData <- function(p, q, K.true, n, loading_means, loading_sd){
	if(missing(p)){p = 40}
	if(missing(q)){q = 4}
	if(missing(K.true)){K.true = 6}
	if(missing(n)){n = 200}
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



	SigmaINV.true <- array(data = 0, dim = c(p,p))
	diag(SigmaINV.true) <- rgamma(p, shape = 1, rate = 1)
	y.true <- array(data = 0, dim = c(n,q))
	Sigma.true <- SigmaINV.true
	diag(Sigma.true) <- 1/diag(SigmaINV.true)
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
		x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true)
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




Kmax <- 20

q.true <- q <- 4
set.seed(10)
p = 40
n = 1000
K = 10
nChains <- 8
dN <- 0.5
#dN <- exp(seq(0,2,length = nChains - 1))
dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax
sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
set.seed(10)
syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )
mod2 = Mclust(syntheticDataset$data, prior = priorControl(functionName="defaultPrior", shrinkage = 0.0), G = 1:20)

sd <- simData(K.true = 10, n = 1000, q = q.true, p = 40)

nChains <- 8
dN <- 0.5
#dN <- exp(seq(0,2,length = nChains - 1))
dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax

outputFolder <- paste0('q.true_',q.true,'q_',q,'_data_',dIter) 
q <- q.true
fabMix( dirPriorAlphas = dirPriorAlphas, rawData = sd$data, outDir = outputFolder, 
        Kmax = Kmax, mCycles = 1200, burnCycles = 200, nIterPerCycle = 10,
	q = q, g = 2, h = 1, alpha_sigma = 2, beta_sigma = 1) 

getStuffForDIC(x_data = originalX, outputFolder = outputFolder, q = q)

dealWithLabelSwitching_same_sigma(x_data = originalX, 
	outputFolder = outputFolder, q = q, 
	compute_regularized_expression = TRUE, Km = Kmax)















