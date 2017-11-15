dir.create("wine")
setwd("wine")
library('fabMix')
library('pgmm')
library('mclust')
data("wine")

########################################################################################################################

#			fabMix 

########################################################################################################################
x<-wine[,-1]
Kmax <- 20      # maximum number of components
nChains <- 8    # number of parallel chains
dN <- 0.5
dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax


myConstant <- 0.5
for(q in 1:5){
	set.seed(1)
	outputFolder <- paste0("fabMix_wine_diff_var_q_", q)
        fabMix( sameSigma = FALSE, dirPriorAlphas = dirPriorAlphas, rawData = x, outDir = outputFolder,
                Kmax = Kmax, mCycles = 1100, burnCycles = 100, q = q, nIterPerCycle = 20,
                g = myConstant, h = myConstant, alpha_sigma = myConstant, beta_sigma = myConstant,  warm_up = 10000) 

	set.seed(1)
	outputFolder <- paste0("fabMix_wine_same_var_q_", q)
        fabMix( sameSigma = TRUE, dirPriorAlphas = dirPriorAlphas, rawData = x, outDir = outputFolder,
                Kmax = Kmax, mCycles = 1100, burnCycles = 100, q = q, nIterPerCycle = 20,
                g = myConstant, h = myConstant, alpha_sigma = myConstant, beta_sigma = myConstant,  warm_up = 1000) 

}

# compute bic values and deal with label switching:
bic <- array(data = NA, dim = c(5,2))
colnames(bic) <- c("different_variance", "same_variance")
for(q in 1:5){
	outputFolder1 <- paste0("fabMix_wine_diff_var_q_", q)
	getStuffForDIC(sameSigma = FALSE, x_data = x, outputFolder = outputFolder1, q = q,  Km = Kmax)
	dealWithLabelSwitching(sameSigma = FALSE, x_data = x, 
		outputFolder =  outputFolder1, q = q, compute_regularized_expression = FALSE, Km = Kmax)

	outputFolder2 <- paste0("fabMix_wine_same_var_q_", q)
	getStuffForDIC(sameSigma = TRUE, x_data = x, outputFolder = outputFolder2, q = q,  Km = Kmax)
	dealWithLabelSwitching(sameSigma = TRUE, x_data = x, 
		outputFolder = outputFolder2, q = q, compute_regularized_expression = FALSE, Km = Kmax)
	bic[q,] <-  c(read.table(paste0(outputFolder1,"/informationCriteria_map_model.txt"))[4,], read.table(paste0(outputFolder2,"/informationCriteria_map_model.txt"))[4,])
}

bestModel <- apply(bic,2,function(y){order(y,decreasing = F)[1]})[order(apply(bic,2,min), decreasing=F)[1]]
# selected parameterization and number of factors:
print(bestModel)
# read single best clustering for this model:
zFABMIX <- as.numeric(read.table("fabMix_wine_diff_var_q_2/singleBestClusterings.txt", header=TRUE)[1,])


########################################################################################################################

#			PGMM 

########################################################################################################################
all_models <- c("CCC", "UUC", "CCU", "CUC", "CUU", "UCC", "UCU", "CCUU", "UCUU", "CUCU","UUCU","UUU")
x <- scale(wine[,-1])
nModels <- length(all_models)
library('foreach')
library('doParallel')

registerDoParallel(cores = nModels)
K1 <- 1
K2 <- 5
KRange <- K1:K2
qRange <- 1:5
myDir <- "pgmm_wine"
dir.create(myDir)
setwd(myDir)
nRandomRuns <- 5
pgmmSeeds <- 100000*(1:nRandomRuns)
foreach(i=1:nModels, .export=ls(envir=globalenv()) ) %dopar% {
		myModel <- all_models[i]
		previousRun <- TRUE
		tryCatch(
	                {
				pgmm_Kmeans <<- pgmmEM(x = x, rG=K1:K2, rq=qRange, icl=FALSE, cccStart=TRUE, zstart=2, loop=1, modelSubset = myModel)
				pgmmRun <<- pgmm_Kmeans
				write.table(pgmmRun$summ_info[[4]], file = paste0("bic_pgmm_",myModel,".txt"))
				write.table(pgmmRun$map, file = paste0("z_pgmm_",myModel,".txt"))
	                        cat("\n\n")
	                }, 
	                error = function(y){ previousRun <<- FALSE; cat(paste0("K-means model ",i,": function failed"),"\n")}
	                )
		for(j in 1:nRandomRuns){
			for(k in KRange){
			for(q in qRange){			
				tryCatch(
					{
						pgmm_Random <<- pgmmEM(x = x, rG=k, rq=q, icl=FALSE, zstart=1, cccStart=TRUE, 
							loop=1, modelSubset = myModel, seed=pgmmSeeds[j])
						if(previousRun == TRUE){
							if( pgmm_Random$summ_info[[4]] > pgmmRun$summ_info[[4]] ){
								pgmmRun <<- pgmm_Random
								write.table(pgmmRun$summ_info[[4]], file = paste0("bic_pgmm_",myModel,".txt"))
								write.table(pgmmRun$map, file = paste0("z_pgmm_",myModel,".txt"))
								write.table(pgmmRun$q, file = paste0("q_pgmm_",myModel,".txt"))
								cat("\n\n")
							}						
						}else{
							previousRun <<- TRUE	
							pgmmRun <<- pgmm_Random
							write.table(pgmmRun$summ_info[[4]], file = paste0("bic_pgmm_",myModel,".txt"))
							write.table(pgmmRun$q, file = paste0("q_pgmm_",myModel,".txt"))
							write.table(pgmmRun$map, file = paste0("z_pgmm_",myModel,".txt"))
							cat("\n\n")
						}
					},
					error = function(y)cat(paste0("random start ",j," for model ",i,": function failed at K = ", k, ", q = ", q),"\n")
					)
			}}
		}

}

# NOTE: the run that gives the best BIC value corresponds to the following command:
#     pgmm_Random <- pgmmEM(x = x, rG=4, rq=4, icl=FALSE, zstart=1, cccStart=TRUE, loop=1, modelSubset = "CUU", seed=100000)
# with printed output:
#     Based on 1 random start, the best model (BIC) for the range of factors and components used is a CUU model with q = 4 and G = 4.
#     The BIC for this model is -11435.43.



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


z_pgmm <- read.table(paste0("z_pgmm_",names(sort(bics, decreasing=TRUE))[1],".txt"))
write.table(z_pgmm, file = "z_pgmm.txt")
setwd("../")

################################################################################################################

# Estimated classifications as reported in Table 4 of the manuscript:

print(table(wine[,1], zFABMIX))
#print(table(wine[,1], EmMfaZ))
print(table(wine[,1], z_pgmm[,1]))













