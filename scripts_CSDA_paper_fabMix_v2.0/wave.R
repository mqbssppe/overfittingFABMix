dir.create("wave")
setwd("wave")


library('fabMix')
library('pgmm')

data("waveDataset1500")
originalX <- waveDataset1500[,-1]

########################################################################################################################

#                       fabMix 

########################################################################################################################
Kmax <- 20              # maximum number of components
nChains <- 8    # number of parallel chains
dN <- 0.25      # NOTE: Larger dN values reduce the acceptance rate of proposed swaps. 
dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax


# The following code will only fit the best model (1 factor and same variance of errors)
set.seed(1)
q = 1
outputFolder <- paste0("wave_same_variance_nFactors_",q)
fabMix( sameSigma = TRUE, dirPriorAlphas = dirPriorAlphas, rawData = originalX, outDir = outputFolder,
        Kmax = Kmax, mCycles = 1600, burnCycles = 100,  q = q,
        g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, warm_up = 10000) 
getStuffForDIC(sameSigma = TRUE, x_data = originalX, outputFolder = outputFolder, q = q, Km = Kmax)
dealWithLabelSwitching(sameSigma = TRUE, x_data = originalX, 
        outputFolder = outputFolder, q = q, compute_regularized_expression = T, Km = Kmax, burn = 100)

zFABMIX <- as.numeric(read.table(paste0(outputFolder,"/singleBestClusterings.txt"), header=TRUE)[1,])


runALL = FALSE
# If you wish to run both parameterizations for a number of factors q = 1,2,3,4 then uncomment the following line:
# runALL = TRUE

if (runALL == TRUE){
for(q in 1:4){
	set.seed(1)
	outputFolder <- paste0("wave_same_variance_nFactors_",q)
	fabMix( sameSigma = TRUE, dirPriorAlphas = dirPriorAlphas, rawData = originalX, outDir = outputFolder,
		Kmax = Kmax, mCycles = 3000, burnCycles = 1000, nIterPerCycle = 10, q = q,
		g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, warm_up = 10000) 
	getStuffForDIC(sameSigma = TRUE, x_data = originalX, outputFolder = outputFolder, q = q, Km = Kmax)
	dealWithLabelSwitching(sameSigma = TRUE, x_data = originalX, 
		outputFolder = outputFolder, q = q, compute_regularized_expression = F, Km = Kmax, burn = 100)
	set.seed(1)
	outputFolder <- paste0("wave_diff_variance_nFactors_",q)
	fabMix( sameSigma = F, dirPriorAlphas = dirPriorAlphas, rawData = originalX, outDir = outputFolder,
		Kmax = Kmax, mCycles = 3000, burnCycles = 1000, nIterPerCycle = 10, q = q,
		g = 0.5, h = 0.5, alpha_sigma = 0.5, beta_sigma = 0.5, warm_up = 10000) 
	getStuffForDIC(sameSigma = F, x_data = originalX, outputFolder = outputFolder, q = q, Km = Kmax)
	dealWithLabelSwitching(sameSigma = F, x_data = originalX, 
		outputFolder = outputFolder, q = q, compute_regularized_expression = F, Km = Kmax, burn = 100)

}
}


########################################################################################################################

#                       PGMM 

########################################################################################################################
x <- scale(originalX, scale=TRUE, center=TRUE)



myDir <- paste0("pgmm")
dir.create(myDir)
setwd(myDir)
library('foreach')
library('fabMix')
library('doParallel')


all_models <- c("CCC", "UUC", "CCU", "CUC", "CUU", "UCC", "UCU", "CCUU", "UCUU", "CUCU","UUCU","UUU")
nModels <- length(all_models)
registerDoParallel(cores = nModels)
K1 <- 1
K2 <- 5	# increase this value to fit greater number of components
KRange <- K1:K2
qRange <- 1:2 # range of the number of factors, adjust as required. (best pgmm model: q = 1, K = 3).

nRuns <- 5 # number of different starts based on random starting values
pgmmSeeds <- (1:nRuns)*100000	# seed per random starting run
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
                pgmmSeeds <- c(100000,200000,300000,400000,500000)
                for(j in 1:5){
                        tryCatch(
                                {
                                        pgmm_Random <<- pgmmEM(x = x, rG=K1:K2, rq=qRange, icl=FALSE, zstart=1, cccStart=TRUE, loop=1, modelSubset = myModel, seed=pgmmSeeds[j])
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
###############################################################################################

# Estimated classifications as reported in Table 3 of the manuscript:

print(table(waveDataset1500$class, zFABMIX))
#print(table(waveDataset1500$class, EmMfaZ))
print(table(waveDataset1500$class, z_pgmm[,1]))


