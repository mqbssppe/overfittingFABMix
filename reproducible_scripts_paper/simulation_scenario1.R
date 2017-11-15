dir.create("simulation_scenario1")
setwd("simulation_scenario1")

library('fabMix')
library('pgmm')
library('mclust')

# Simulation set-up:

n = 500               # sample size
p = 40                # number of variables
q = 4                 # number of factors
K = 10                # number of clusters


#	simulate data as described in the Appendix A.2 for scenario 1 (default setting)
sINV_diag = 1/((1:p))	 				# diagonal of inverse variance of errors
set.seed(10)
syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag)
qRange <- 1:8	# range of values for the number of factors
########################################################################################################################

#                       fabMix 

########################################################################################################################

Kmax <- 20		# number of components for the overfitted mixture model
nChains <- 8		# number of parallel heated chains
dN <- 2         	# parameter the controls difference between successive chains
dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax	# Dirichlet prior parameter per chain


for(q in qRange){
	set.seed(1)
		outputFolder <- paste0("fabMix_scenario1_sameSigma_q_",q)
		fabMix(sameSigma = TRUE, dirPriorAlphas = dirPriorAlphas, rawData = syntheticDataset$data, outDir = outputFolder,
		        Kmax = Kmax, mCycles = 1600, burnCycles = 100, q = q,
		        g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
		        warm_up_overfitting = 500, warm_up = 4000) 
		getStuffForDIC(sameSigma = TRUE, x_data = syntheticDataset$data, outputFolder = outputFolder, q = q, Km = Kmax)
		dealWithLabelSwitching(sameSigma = TRUE, x_data = syntheticDataset$data, 
		        outputFolder = outputFolder, q = q, compute_regularized_expression = FALSE, Km = Kmax)
	set.seed(1)
		outputFolder <- paste0("fabMix_scenario1_diffSigma_q_",q)
		fabMix(sameSigma = FALSE, dirPriorAlphas = dirPriorAlphas, rawData = syntheticDataset$data, outDir = outputFolder,
		        Kmax = Kmax, mCycles = 1600, burnCycles = 100, q = q,
		        g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, 
		        warm_up = 1000) 
		getStuffForDIC(sameSigma = FALSE, x_data = syntheticDataset$data, outputFolder = outputFolder, q = q, Km = Kmax)
		dealWithLabelSwitching(sameSigma = FALSE, x_data = syntheticDataset$data, 
		        outputFolder = outputFolder, q = q, compute_regularized_expression = FALSE, Km = Kmax)
}

# read bic values:
bic <- array(data = NA, dim = c(length(qRange),2))
colnames(bic) <- c("sameSigma", "diffSigma")
for(q in qRange){
        outputFolder1 <- paste0("fabMix_scenario1_sameSigma_q_", q)
        outputFolder2 <- paste0("fabMix_scenario1_diffSigma_q_", q)
        bic[q,] <-  c(read.table(paste0(outputFolder1,"/informationCriteria_map_model.txt"))[4,], read.table(paste0(outputFolder2,"/informationCriteria_map_model.txt"))[4,])
}

bestModel <- apply(bic,2,function(y){order(y,decreasing = F)[1]})[order(apply(bic,2,min), decreasing=F)[1]]
# selected parameterization and number of factors (it corresponds to q=4 factors with same variance of errors per cluster):
print(bestModel)
# read single best clustering for this model:
model_parameterization <- names(bestModel)
nFactors <- as.numeric(bestModel)
cat(paste0("Selected parameterization: ", model_parameterization, ", selected number of factors: ", nFactors), "\n")
zFABMIX <- as.numeric(read.table(paste0("fabMix_scenario1_",model_parameterization,"_q_",nFactors,"/singleBestClusterings.txt"), header=TRUE)[1,])


########################################################################################################################

#                       PGMM 

########################################################################################################################

library('foreach')
library('fabMix')
library('doParallel')

x <- scale(syntheticDataset$data, scale=TRUE, center=TRUE)

all_models <- c("CCC", "UUC", "CCU", "CUC", "CUU", "UCC", "UCU", "CCUU", "UCUU", "CUCU","UUCU","UUU")
nModels <- length(all_models)
registerDoParallel(cores = nModels)
K1 <- 1
K2 <- K + 5	# max number of components: in this example K2 = 15. Recall that K = 10 denotes the true number of clusters. 
KRange <- K1:K2
myDir <- "pgmm_scenario1"
dir.create(myDir)
setwd(myDir)
nRandomRuns <- 5
pgmmSeeds <- (1:nRandomRuns)*100000
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
                for(j in 1:nRandomRuns){
                        tryCatch(
                                {
                                        pgmm_Random <<- pgmmEM(x = x, rG=K1:K2, rq=qRange, icl=FALSE, zstart=1, cccStart=TRUE, loop=1, modelSubset = myModel, seed=pgmmSeeds[j])
                                        if(previousRun == TRUE){
                                                if( pgmm_Random$summ_info[[4]] > pgmmRun$summ_info[[4]] ){
                                                        pgmmRun <<- pgmm_Random
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
# read best bics per model
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

# best pgmm clustering
z_pgmm <- read.table(paste0("z_pgmm_",names(sort(bics, decreasing=TRUE))[1],".txt"))
write.table(z_pgmm, file = "z_pgmm.txt")
setwd("../")

################################################################################################################

# Estimated AdjustedRandIndex and Number of clusters as reported in Table 2-scenario1 of the manuscript:

scenario1_result <- array(data = NA, dim = c(3,2))
colnames(scenario1_result) <- c("estimated_number_of_clusters", "adjusted_rand_index")
rownames(scenario1_result) <- c("fabMix", "EMMIXmfa", "pgmm")

scenario1_result[1, ] <- c(length(table(zFABMIX)), adjustedRandIndex(syntheticDataset$class, zFABMIX))
#scenario1_result[2, ] <- c(length(table(EmMfaZ)), adjustedRandIndex(syntheticDataset$class, EmMfaZ))
scenario1_result[3, ] <- c(length(table(z_pgmm[,1])), adjustedRandIndex(syntheticDataset$class, z_pgmm[,1]))

print(scenario1_result)






