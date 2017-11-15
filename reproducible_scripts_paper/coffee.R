dir.create("coffee")
setwd("coffee")

library("pgmm")
library('fabMix')
library('mclust')
data("coffee")
x<-coffee[,-c(1,2)]
x<-scale(x)

########################################################################################################################

#                       fabMix 

########################################################################################################################

Kmax <- 20              # maximum number of components
nChains <- 8    # number of parallel chains
dN <- 1         
dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax

c1 <- 0.5
c2 <- 0.5

set.seed(1)
for (q in 1:5){
	outputFolder <- paste0("diff_var_q_", q, "_c1_", c1, "_c2_", c2)
	fabMix( sameSigma = FALSE, dirPriorAlphas = dirPriorAlphas, rawData = coffee[,-c(1,2)], outDir = outputFolder,
		Kmax = Kmax, mCycles = 2500, burnCycles = 500,  q = q,
		g = c1, h = c2, alpha_sigma = c1, beta_sigma = c2) 
	getStuffForDIC(sameSigma = FALSE, x_data = coffee[,-c(1,2)], outputFolder = outputFolder, q = q, Km = Kmax)
	dealWithLabelSwitching(sameSigma = FALSE, x_data = coffee[,-c(1,2)], 
		outputFolder = outputFolder, q = q, compute_regularized_expression = FALSE, Km = Kmax)
}


set.seed(1)
for (q in 1:5){
	outputFolder <- paste0("same_var_q_", q, "_c1_", c1, "_c2_", c2)
	fabMix( sameSigma = TRUE, dirPriorAlphas = dirPriorAlphas, rawData = coffee[,-c(1,2)], outDir = outputFolder,
		Kmax = Kmax, mCycles = 2500, burnCycles = 500,  q = q,
		g = c1, h = c2, alpha_sigma = c1, beta_sigma = c2) 
	getStuffForDIC(sameSigma = TRUE, x_data = coffee[,-c(1,2)], outputFolder = outputFolder, q = q, Km = Kmax)
	dealWithLabelSwitching(sameSigma = TRUE, x_data = coffee[,-c(1,2)], 
		outputFolder = outputFolder, q = q, compute_regularized_expression = FALSE, Km = Kmax)
}



bic <- array(data = NA, dim = c(5,2))
colnames(bic) <- c("different_variance", "same_variance")

bic[1,1] <- read.table("diff_var_q_1_c1_0.5_c2_0.5/informationCriteria_map_model.txt")[4,]
bic[2,1] <- read.table("diff_var_q_2_c1_0.5_c2_0.5/informationCriteria_map_model.txt")[4,]
bic[3,1] <- read.table("diff_var_q_3_c1_0.5_c2_0.5/informationCriteria_map_model.txt")[4,]
bic[4,1] <- read.table("diff_var_q_4_c1_0.5_c2_0.5/informationCriteria_map_model.txt")[4,]
bic[5,1] <- read.table("diff_var_q_5_c1_0.5_c2_0.5/informationCriteria_map_model.txt")[4,]

bic[1,2] <- read.table("same_var_q_1_c1_0.5_c2_0.5/informationCriteria_map_model.txt")[4,]
bic[2,2] <- read.table("same_var_q_2_c1_0.5_c2_0.5/informationCriteria_map_model.txt")[4,]
bic[3,2] <- read.table("same_var_q_3_c1_0.5_c2_0.5/informationCriteria_map_model.txt")[4,]
bic[4,2] <- read.table("same_var_q_4_c1_0.5_c2_0.5/informationCriteria_map_model.txt")[4,]
bic[5,2] <- read.table("same_var_q_5_c1_0.5_c2_0.5/informationCriteria_map_model.txt")[4,]

bestModel <- apply(bic,2,function(y){order(y,decreasing = F)[1]})[order(apply(bic,2,min), decreasing=F)[1]]
# selected parameterization and number of factors:
print(bestModel)
# read single best clustering for this model:
zFABMIX <- as.numeric(read.table("diff_var_q_1_c1_0.5_c2_0.5/singleBestClusterings.txt", header=TRUE)[1,])

########################################################################################################################

#                       PGMM 

########################################################################################################################

x<-coffee[,-c(1,2)]
x<-scale(x)

myDir <- paste0("pgmm_1000")
dir.create(myDir)
setwd(myDir)
library(pgmm)
library(foreach)
library(fabMix)
library(doParallel)


all_models <- c("CCC", "UUC", "CCU", "CUC", "CUU", "UCC", "UCU", "CCUU", "UCUU", "CUCU","UUCU","UUU")
nModels <- length(all_models)
registerDoParallel(cores = nModels)
K1 <- 1
K2 <- 5
KRange <- K1:K2
qRange <- 1:5
nRuns <- 1000
foreach(i=1:nModels, .export=ls(envir=globalenv()) ) %dopar% {
		myModel <- all_models[i]
		previousRun <- TRUE
	        tryCatch(
                        {
				pgmm_Kmeans <<- pgmmEM(x = x, rG=K1:K2, rq=qRange, icl=FALSE, cccStart=TRUE, zstart=2, loop=1, modelSubset = myModel)
				pgmmRun <<- pgmm_Kmeans
				write.table(pgmmRun$summ_info[[4]], file = paste0("bic_pgmm_",myModel,".txt"))
				write.table(pgmmRun$map, file = paste0("z_pgmm_",myModel,".txt"))
				write.table(pgmmRun$q, file = paste0("q_",myModel,".txt"))
                                cat("\n\n")
                        }, 
                        error = function(y){ previousRun <<- FALSE; cat(paste0("K-means model ",i,": function failed"),"\n")}
                        )
		pgmmSeeds <- (1:nRuns)*100000
		for(j in 1:nRuns){
	        	tryCatch(
				{
					pgmm_Random <<- pgmmEM(x = x, rG=K1:K2, rq=qRange, icl=FALSE, zstart=1, cccStart=TRUE, loop=1, modelSubset = myModel, seed=pgmmSeeds[j])
					cat('\n')
					if(previousRun == TRUE){
						if( pgmm_Random$summ_info[[4]] > pgmmRun$summ_info[[4]] ){
							pgmmRun <<- pgmm_Random
							cat('\n')
							       cat('***********************************************************************************','\n')							
							cat(paste0('*** found better estimates for model parameterization: ', myModel, ', with: q = ',pgmmRun$q , ' and K = ', pgmmRun$g), '\n')
							       cat('***********************************************************************************','\n')							
							write.table(pgmmRun$summ_info[[4]], file = paste0("bic_pgmm_",myModel,".txt"))
							write.table(pgmmRun$map, file = paste0("z_pgmm_",myModel,".txt"))
							write.table(pgmmRun$q, file = paste0("q_",myModel,".txt"))
						        cat("\n\n")
						}						
					}else{
						previousRun <<- TRUE	
						pgmmRun <<- pgmm_Random
						write.table(pgmmRun$summ_info[[4]], file = paste0("bic_pgmm_",myModel,".txt"))
						write.table(pgmmRun$map, file = paste0("z_pgmm_",myModel,".txt"))
						write.table(pgmmRun$q, file = paste0("q_",myModel,".txt"))
				                cat("\n\n")
					}
				},
				error = function(y)cat(paste0("random start ",j," for model ",i,": function failed"),"\n")
				)
		}

}


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

z_pgmm <- read.table(paste0("z_pgmm_",names(sort(bics, decreasing=TRUE))[1],".txt"))[,1]
write.table(z_pgmm, file = "z_pgmm.txt")
setwd("../")



################################################################################################################

# Estimated classifications as reported in Table 5 of the manuscript:

print(table(coffee[,1], zFABMIX))
#table(coffee[,1], EmMfaZ)
print(table(coffee[,1], z_pgmm))
















