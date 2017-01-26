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

Kmax <- 20
q <- q.true
#qIndex <- (1:10)
#for(q in qIndex){

	nChains <- 8
        dN <- 2.5
	dN <- exp(seq(0,2,length = nChains - 1))
        dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax

	outputFolder <- paste0('q.true_',q.true,'q_',q,'_data_',dIter) 

	fabMix( dirPriorAlphas = dirPriorAlphas, rawData = originalX, outDir = outputFolder, 
	        Kmax = Kmax, mCycles = 1200, burnCycles = 200, nIterPerCycle = 20,
		q = q, g = 2, h = 1, alpha_sigma = 2, beta_sigma = 1) 

	getStuffForDIC(x_data = originalX, outputFolder = outputFolder, q = q)

	dealWithLabelSwitching_same_sigma(x_data = originalX, 
		outputFolder = outputFolder, q = q, 
		compute_regularized_expression = TRUE, Km = Kmax)

#}




