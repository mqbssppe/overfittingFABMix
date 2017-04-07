# *******************************************************************************
# (a) Simulate a dataset along the lines of the paper.
#       It follows the scenario used in subsection:  
#       ``Estimating the number of clusters''                   
# (b) Run fabMix + EMMIXmfa + pgmm (mixtures of factor analyzers)
# (c) Compare estimated clusters against the true classification 
# *******************************************************************************


library('fabMix')
library('pgmm')
library('EMMIXmfa')
library('mclust')	# for the adjusted rand index

#_______________________________________________________________________________________________________________________
#__________________________(a): Data simulation_________________________________________________________________________

# Simulation parameters							

n = 500                 # sample size ***NOTE: in the paper we have used n = 1000 and n = 2000***
p = 40                  # number of variables
q = 4                   # number of factors
K = 10                  # number of clusters
sINV_diag = 1/((1:p))   # diagonal of inverse variance of errors
set.seed(10)
syntheticDataset <- simData(K.true = K, n = n, q = q, p = p, sINV_values = sINV_diag )

#________________________________________________________________________________________________________________________
#__________________________(b): Fit MFA models___________________________________________________________________________

Kmax <- 20		# maximum number of components

# All MFA methods will run assuming that the number of factors is equal to the true value (q = 4).
# For a complete analysis, you should consider various values of q (eg: q = 1, 2,... , 10) 
#	and select the model with the smallest BIC or AIC. 

#	(b.1): pgmm package .............................................................................................
#	Run pgmm for a range of MFA models with K = 1,..., Kmax and select the best model with BIC			#
set.seed(1)														#
pgmmRun <- pgmmEM(x = syntheticDataset$data, rG=1:Kmax, rq=4, class=NULL, icl=FALSE, zstart=2, 				#
		cccStart=TRUE, loop=1, zlist=NULL, modelSubset=c("UCU"))						#
K.pgmm <- pgmmRun$g													#
ARI.pgmm <- adjustedRandIndex(syntheticDataset$class, apply(pgmmRun$zhat,1,function(y){order(y,decreasing=T)[1]}))	#
#.......................................................................................................................#
#	(b.2): EMMIXmfa package..........................................................................................
#	Run EMMIXmfa for a range of MFA models with K = 1,..., Kmax and select the best model with BIC			#
set.seed(1)														#
EmMfaRun4 <- vector("list", length = Kmax)										#
for(k in 1:Kmax){													#
        EmMfaRun4[[k]] <- mfa(Y = syntheticDataset$data, g=k, q=4, itmax=1000, nkmeans=1, nrandom=5, 			#
	sigmaType = "unique", Dtype = "common")										#
}															#
K.EMMIXmfa <- order(unlist(lapply(EmMfaRun4, function(y)y$BIC)), decreasing=F)[1]					#
ARI.EMMIXmfa <- adjustedRandIndex(EmMfaRun4[[K.EMMIXmfa]]$clust, syntheticDataset$class)				#
#.......................................................................................................................#
#	(b.3): fabMix package............................................................................................
# 	fabMix parameters: 												#
nChains <- 8  	# number of parallel chains										#
dN <- 3		# NOTE: Larger dN values reduce the acceptance rate of proposed swaps. 					#
		#	The user should carefully examine this parameter to achieve 					#
		#	reasonable acceptance rates. In the paper values 1 <= dN <= 5 are considered. 			#
# Dirichlet prior of mixture weights per chain.										#
#   The target chain corresponds to the first entry.									#
dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/Kmax									#
outputFolder <- paste0("fabMixExample_nFactors_",q)									#
# step.1) Run algorithm (warning: 8 parallel threads needed)								#
fabMix( dirPriorAlphas = dirPriorAlphas, rawData = syntheticDataset$data, outDir = outputFolder, 			#
	Kmax = Kmax, mCycles = 1200, burnCycles = 200, q = q) 								#
# step.2) Compute information criteria for the given value of q								#
getStuffForDIC(x_data = syntheticDataset$data, outputFolder = outputFolder, q = q)					#
# step.3) Undo the label switching for the most probable number of alive clusters					#
dealWithLabelSwitching_same_sigma(x_data = syntheticDataset$data, 							#
	outputFolder = outputFolder, q = q, compute_regularized_expression = FALSE, Km = Kmax)				#
#	Get estimated clusters using the ECR reordering algorithm:							#
fabMix_clusters <- as.matrix(read.table(paste0(outputFolder,"/singleBestClusterings.txt"), header= TRUE))[1,]		#
ARI.fabMix <- adjustedRandIndex(syntheticDataset$class, fabMix_clusters)						#
K.fabMix <- length(table(fabMix_clusters))										#
#.......................................................................................................................#


# compute adjusted rand index with respect to the true class:
adjRandIndex <- array(data = NA, dim = c(3,2))
rownames(adjRandIndex) <- c("EMMIXmfa", "fabMix", "pgmm")
colnames(adjRandIndex) <- c("ARI", "nClusters")
adjRandIndex["fabMix",1] <-  ARI.fabMix
adjRandIndex["EMMIXmfa",1] <-  ARI.EMMIXmfa
adjRandIndex["pgmm",1] <-  ARI.pgmm
adjRandIndex["fabMix",2] <-  K.fabMix
adjRandIndex["EMMIXmfa",2] <-  K.EMMIXmfa
adjRandIndex["pgmm",2] <- K.pgmm

# the result should be similar to the following output:
> adjRandIndex
               ARI nClusters
EMMIXmfa 0.8432545        11
fabMix   1.0000000        10
pgmm     0.7990348         9

