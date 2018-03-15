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
x <- scale(wine[,-1])
myDir <- "pgmm_wine"
dir.create(myDir)
setwd(myDir)

# 1 run based on k-means starting values:
pgmm_Kmeans <- pgmmEM(x = x, rG=1:5, rq=1:5, relax = TRUE, icl=FALSE, cccStart=TRUE, zstart=2, loop=1)

# 5 runs based on random starting values:
wine_pgmmEM = pgmmEM(x,rG=1:5,rq=1:5,relax=TRUE,zstart=1,loop=5)
# best BIC value corresponds to wine_pgmmEM
z_pgmm <- wine_pgmmEM$map
setwd("../")

################################################################################################################

# Estimated classifications as reported in Table 4 of the manuscript:

print(table(wine[,1], zFABMIX))
#print(table(wine[,1], EmMfaZ))
print(table(wine[,1], z_pgmm))













