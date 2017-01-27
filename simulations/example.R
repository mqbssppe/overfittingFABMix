#The following data generation mechanism is used:

library(fabMix)
library(flexmix)
library(mclust)

n = 1000              # sample size
p = 40                # number of variables
q = 4                 # number of factors
K = 10                # number of clusters
#			scenario 1:
sINV_diag = 1/((1:p)) # diagonal of inverse variance of errors
#			or you can also use scenario 2:
#sINV_diag = rgamma(p, shape = 1, rate = 1)
#			note that scenario is more challenging
#			as it produces more noisy observations
#			for features j --> p 
set.seed(10)
syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )
#Simulation parameters: 
   n = 1000   # sample size
   p = 40     # number of variables
   q = 4      # number of factors
   K = 10     # number of  clusters
Kmax <- 20    # number of overfitted mixture components
nChains <- 8  # number of parallel chains
dN <- 1
# Dirichlet prior of mixture weights per chain.
#   The target chain corresponds to the first entry.
dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax


outputFolder <- paste0("fabMixExample_nFactors_",q)
# step.1) Run algorithm (warning: 8 parallel threads needed, computation takes ~ 1.5 - 2 hours)
ptm <- proc.time()
fabMix( dirPriorAlphas = dirPriorAlphas, 
        rawData = syntheticDataset$data, 
        outDir = outputFolder, Kmax = Kmax, mCycles = 1200, 
        burnCycles = 200, q = q) 
timeNeeded <- proc.time() - ptm
# step.2) Compute information criteria for the given value of q
getStuffForDIC(x_data = syntheticDataset$data, outputFolder = outputFolder, q = q)

# step.3) Undo the label switching for the most probable number of alive clusters
dealWithLabelSwitching_same_sigma(x_data = syntheticDataset$data, 
	outputFolder = outputFolder, q = q, 
	compute_regularized_expression = TRUE, Km = Kmax)

# for a complete analysis, you should rerun steps 1, 2, 3 with various values of q (eg: q = 1, 2,... , 10) 
#	and select the model with the smallest BIC or AIC. It corresponds to q = 4. 


fabMix_clusters <- as.matrix(read.table(paste0(outputFolder,"/singleBestClusterings.txt"), header= TRUE))[1,]


# Run mclust

mclustModel = Mclust(syntheticDataset$data, prior = priorControl(functionName="defaultPrior"), G = 1:20)
mclust_clusters <- mclustModel$classification

# Run flexmix
flexmixModel <- initFlexmix(syntheticDataset$data ~ 1, k = 1:20, model = FLXMCmvnorm(diagonal = FALSE), nrep = 10)
getFlexmixResult <- getModel(flexmixModel, which = "BIC")
flexmix_clusters <- getFlexmixResult@cluster

# compute adjusted rand index with respect to the true class:
adjRandIndex <- array(data = NA, dim = c(3,2))
rownames(adjRandIndex) <- c("fabMix", "flexmix", "mclust")
colnames(adjRandIndex) <- c("ARI", "nClusters")
adjRandIndex["fabMix",1] <-  adjustedRandIndex(syntheticDataset$class, fabMix_clusters)
adjRandIndex["flexmix",1] <-  adjustedRandIndex(syntheticDataset$class, flexmix_clusters)
adjRandIndex["mclust",1] <-  adjustedRandIndex(syntheticDataset$class, mclust_clusters)
adjRandIndex["fabMix",2] <-  length(table(fabMix_clusters))
adjRandIndex["flexmix",2] <-  length(table(flexmix_clusters))
adjRandIndex["mclust",2] <-  length(table(mclust_clusters))

# the result should be similar to the following output:

> adjRandIndex
              ARI nClusters
fabMix  1.0000000        10
flexmix 0.3600983         3
mclust  0.7068138        19









