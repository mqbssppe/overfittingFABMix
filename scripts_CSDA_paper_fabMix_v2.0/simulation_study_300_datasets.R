# This script will take toooooo long to complete (weeks). 
# If possible, split the script to independent jobs and run each one on a server. 
# This script demands 8 cores for fabMix, and 12 cores for pgmm. 

library('fabMix')
library('pgmm')
library('mclust') #   for computing the adjustedRandIndex
library('foreach')
library('doParallel')


dir.create("simulation_study_300_datasets")
setwd("simulation_study_300_datasets")

p = 40                # number of variables
q.true = 4            # true number of factors


# default values used for the analysis in the paper:
mySeeds <- c(10, 30303, 38484, 555, 2020, 5783, 99999, 4348441, 1234, 435969)
nSeeds <- length(mySeeds)	# we generated 10 datasets per case, using the seeds defined in mySeeds
nRange <- c(500, 1000, 2000)	# sample size values
K.trueRange <- 1:10		# true number of clusters values
mCycles <- 1600		# mcmc cycles for fabMix after warm_up iterations. 1 cycle = 10 iterations
warm_up <- 4000		# iterations that correspond to warm up period 
qRange <- q.true	# Range for the number of factors.
			# 	NOTE: only the model with the true number of factors is fitted.
			# 	If you want to consider various values for the number of factors, uncomment the following line:
# qRange <- 1:(2*q.true)


# results are stored in the following data frame:
results <- data.frame( K.true = numeric(),		# true number of clusters
			n = numeric(),			# sample size
			method = character(),		# name of each method
			ARI = numeric(),		# adjusted Rand Index
			K.estimate = numeric()		# estimated number of clusters
		)
# The `results` data frame is also saved on the file 'results.txt' while the loop runs. 

for(dIter in 1:nSeeds){			# dataset: 1,2, ..., 10
	for(n in nRange){		# sample size: 500, 1000, 2000
		for(K in K.trueRange){	# true number of clusters: 1,2,...,10

#			Simulate dataset according to scenario 1, with K = K.true, q = q.true 
			set.seed(mySeeds[dIter])
			sINV_diag = 1/((1:p))	# diagonal of the inverse matrix of variance of errors
			syntheticDataset <- simData(K.true = K, n = n, q = q.true, p = p, sINV_values = sINV_diag )
#			FABMIX 
			Kmax <- 20      # number of components for the overfitted MFA
			nChains <- 8	# number of heated chains
			dN <- 1         # temperature parameter
			dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax	# temperatures
			bics <- numeric(length(qRange))
			names(bics) <- qRange
			for(q in qRange){
					outputFolder <- paste0("fabMix_sameSigma_n",n,"_K",K,"_data",dIter,"_q_",q)
					fabMix(sameSigma = TRUE, dirPriorAlphas = dirPriorAlphas, rawData = syntheticDataset$data, outDir = outputFolder,
						Kmax = Kmax, mCycles = mCycles, burnCycles = 100, q = q,
						g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
						warm_up = warm_up) 
					getStuffForDIC(sameSigma = TRUE, x_data = syntheticDataset$data, outputFolder = outputFolder, q = q)
					dealWithLabelSwitching(sameSigma = TRUE, x_data = syntheticDataset$data, 
						outputFolder = outputFolder, q = q, compute_regularized_expression = FALSE, Km = Kmax)
					bics[as.character(q)] <- read.table(paste0(outputFolder, "/informationCriteria_map_model.txt"))[4,]
			}
			best_q <- names(bics)[order(bics, decreasing = F)[1]]	# selected number of factors
			fileExists <- file.exists(paste0("fabMix_sameSigma_n",n,"_K",K,"_data",dIter,"_q_",best_q,"/singleBestClusterings.txt"))
			if(fileExists){
				# case with K_estimated > 1
				z_FABMIX <- as.numeric(read.table(paste0("fabMix_sameSigma_n",n,"_K",K,"_data",dIter,"_q_",best_q,"/singleBestClusterings.txt"),header=TRUE)[1,])	# estimated clusters according to ECR algorithm
			}else{
				# case with K_estimated = 1 (our algorithm does not produce this file in this case, so I set it manually here:
				z_FABMIX <- rep(1,n)
			}
			ari <- adjustedRandIndex(syntheticDataset$class, z_FABMIX)
			kEstimate <- length(table(z_FABMIX))
			tmp <- data.frame(K.true = K, n = n, method="fabMix", ARI = ari, K.estimate = kEstimate)
			results <- rbind(results, tmp)
#			FABMIX end
#			PGMM 
			x <- scale(syntheticDataset$data, scale=TRUE, center=TRUE)
			all_models <- c("CCC", "UUC", "CCU", "CUC", "CUU", "UCC", "UCU", "CCUU", "UCUU", "CUCU","UUCU","UUU")
			nModels <- length(all_models)
			registerDoParallel(cores = nModels)
			K1 <- max(c(1, K-2))
			K2 <- K+5
			KRange <- K1:K2
			myDir <- paste0("n_",n,"K_",K,"_data_",dIter)
			dir.create(myDir)
			setwd(myDir)
			nRandomRuns <- 5
			pgmmSeeds <- 100000*(1:nRandomRuns)
			foreach(i=1:nModels, .export=ls(envir=globalenv()) ) %dopar% {
					myModel <- all_models[i]
					previousRun <- TRUE
					tryCatch(
						{
							pgmm_Kmeans <<- pgmmEM(x = x, rG=K1:K2, rq=q, icl=FALSE, cccStart=TRUE, zstart=2, loop=1, modelSubset = myModel)
							pgmmRun <<- pgmm_Kmeans
							write.table(pgmmRun$summ_info[[4]], file = paste0("bic_pgmm_",myModel,".txt"))
							write.table(pgmmRun$map, file = paste0("z_pgmm_",myModel,".txt"))
							cat("\n\n")
						}, 
						error = function(y){ previousRun <<- FALSE; cat(paste0("K-means model ",i,": function failed"),"\n")}
						)
					for(j in 1:nRandomRuns){
						tryCatch(
							{
								pgmm_Random <<- pgmmEM(x = x, rG=K1:K2, rq=q, icl=FALSE, zstart=1, cccStart=TRUE, loop=1, modelSubset = myModel, seed=pgmmSeeds[j])
								if(previousRun == TRUE){
								        if( pgmm_Random$summ_info[[4]] > pgmmRun$summ_info[[4]] ){
								                pgmmRun <<- pgmm_Random
								                write.table(pgmmRun$summ_info[[4]], file = paste0("bic_pgmm_",myModel,".txt"))
								                write.table(pgmmRun$map, file = paste0("z_pgmm_",myModel,".txt"))
								                cat("\n\n")
								        }                                               
								}else{
								        previousRun <<- TRUE    
								        pgmmRun <<- pgmm_Random
								        write.table(pgmmRun$summ_info[[4]], file = paste0("bic_pgmm_",myModel,".txt"))
								        write.table(pgmmRun$map, file = paste0("z_pgmm_",myModel,".txt"))
								        cat("\n\n")
								}
							},
							error = function(y)cat(paste0("random start ",j," for model ",i,": function failed"),"\n")
							)
					}

			}
			stopImplicitCluster()
			bics <- rep(NA, nModels)
			names(bics) <-  all_models
			i <- 0
			for (model in all_models){
				i <- i + 1
				myFile <- paste0("bic_pgmm_",model,".txt")

				if(file.exists(myFile)){
				        bics[i] <- read.table(myFile)[1,1]
				}       
			}
			myFile <-paste0("z_pgmm_",names(sort(bics, decreasing=TRUE))[1],".txt")	
			if(file.exists(myFile)){
				z_pgmm <- read.table(myFile)[,1]
				write.table(z_pgmm, file = "z_pgmm.txt")
			}else{z_pgmm <- rep(NA, n)}
			setwd("../")
			ari <- adjustedRandIndex(syntheticDataset$class, z_pgmm)
			kEstimate <- length(table(z_pgmm))
			tmp <- data.frame(K.true = K, n = n, method="pgmm", ARI = ari, K.estimate = kEstimate)
			results <- rbind(results, tmp)
#			PGMM end
			write.table(results, file = "results.txt")	# saves the results while the loop runs...
		}	
	}
}


# visualize results:
# number of clusters
par(mfrow = c(2,length(nRange)))
for (n in nRange){
	ind <- which(results$n == n)
	boxplot(results[ind, ]$K.estimate ~ interaction( results[ind, ]$method, results[ind, ]$K.true ), border = c("red", "green"), ylab = "number of clusters", xlab = "method", main = paste0("n = ", n))
	points(1:(length(K.trueRange)*2), rep(K.trueRange, each = 2), col = 1, pch = 3, cex= 2)
	legend("topleft", col = c("red", "green", "black"), c("fabMix", "pgmm", "K true"), pch = c(0,0,3), lty = c(1,1,0))
}
# adjusted rand index
for (n in nRange){
	ind <- which(results$n == n)
	boxplot(results[ind, ]$ARI ~ interaction( results[ind, ]$method, results[ind, ]$K.true ), border = c("red", "green"), ylab = "adjusted rand index", xlab = "method", main = paste0("n = ", n))
	legend("bottomleft", col = c("red", "green"), c("fabMix", "pgmm"), pch = c(0,0), lty = c(1,1))
}



