source('~/Dropbox/sparseFA_MIX/heated_prior/q_0_same_sigma.R')
	q <- 0
	heated_chains( rawData = originalX, outDir = paste0('q_',q),Kmax = Kmax, mCycles = 1100, burnCycles = 100, g = 2, h = 0.001, alpha_sigma = 2, beta_sigma = 0.001, q = q)
	setwd("q_0/")
	library('label.switching')
	z <- as.matrix(read.table("zValues.txt"))
	logl <- read.table("kValues.txt", header=T)
	burn <- 1:1
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
	ls <- label.switching(method = "ECR", zpivot = zPivot, z = z, K = K)
	z_Est <- ls$clusters[1,]
	write.table(file = "z_Est.txt", z_Est, row.names = FALSE, col.names=FALSE)
setwd("../")

source('~/Dropbox/sparseFA_MIX/heated_prior/csf_tdp_same_sigma.R')
for(q in 1:10){
	p <- dim(originalX)[2]
	v_r <- numeric(p) 
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}
	outDir <- paste0('q_',q)
	heated_chains( rawData = originalX, outDir = outDir, mCycles = 1500, burnCycles = 500, g = 2, h = 0.001, alpha_sigma = 2, beta_sigma = 0.001, q = q, zStart = z_Est)  
	setwd("outDir")
		z <- as.matrix(read.table("zValues.txt"))
		logl <- read.table("kValues.txt", header=T)
		burn <- 1:1
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
		ls <- label.switching(method = "ECR", zpivot = zPivot, z = z, K = K)
		z_Est <- ls$clusters[1,]
		write.table(file = "z_Est.txt", z_Est, row.names = FALSE, col.names=FALSE)
	setwd("../")
}

