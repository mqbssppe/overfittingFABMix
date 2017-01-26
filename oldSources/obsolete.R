
heated_chains_omega <- function( alpha_prior, rawData, outDir, Kmax, mCycles, burnCycles, g_h_parameters, alpha_sigma, beta_sigma, q, normalize, thinning, zStart){
	if(missing(Kmax)){Kmax <- 20}
	if(missing(zStart)){zStart = FALSE}
	if( missing( alpha_prior) ){
		alpha_prior <- rep( 1/Kmax, Kmax)
	}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g_h_parameters)){
		g_h_parameters <- array(data = NA, dim = c(4,2))
		h_seq <- 10^(-seq(3,4, length = 4))	# auto + (99) leitourgei OK alla acceptance rate ~ 1.
#		h_seq <- 10^(-seq(3,6, length = 8))	# + step (99)
#		g_seq <- 2 + 2*(0:8)	# + step (66) alla einai gia ton poutso
		for(i in 1:4){
 			g_h_parameters[i, ] <- c(2, h_seq[i])   ######  (99)
#			g_h_parameters[i, ] <- c(g_seq[i], 0.001) ##### (66)
		}
	}
	nChains <- dim(g_h_parameters)[1]
	if(missing(alpha_sigma)){alpha_sigma <- 0.1}
	if(missing(beta_sigma)){beta_sigma <- 0.1}

	dir.create(outDir)
	setwd(outDir)
	registerDoParallel(cores = nChains)
	outputDirs <- paste0('alpha_',1:nChains)
	originalX <- rawData
	x_data <- originalX
	if( missing(thinning) ){thinning = 10}
	if( thinning < 1 ){ stop('thinning should be larger than or equal to 1.') }
	thinning <- floor(thinning)
	if( missing(normalize) ){normalize <- TRUE}
	cat('\n')
	cat(paste0("-    p = ", p, ", q = ", q, ", n = ",n,", g = ", g_h_parameters[1,1], ", h = ", g_h_parameters[1,2], ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	cat(paste0('-    Using Nchains = ', nChains),'\n')
	if( normalize == TRUE ){
		x_data <- scale(originalX, center = TRUE, scale = TRUE)
		cat('-    The sampler uses standardized data.','\n')
	}
	if( normalize == FALSE ){
		x_data <- rawData
		cat('-    The sampler uses raw data.','\n')
	}
	kValues <- array(data = NA, dim = c(mCycles, nChains))
	mh_acceptance_rate <- 0
	dir.create('tmpDir')
	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}


	#	initialization
	iteration <- 1
	foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		truncatedDirichletProcess(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
			Kmax = Kmax, m = 5, thinning = 1, burn = 4, alpha_prior= alpha_prior, g = g_h_parameters[myChain,1], h = g_h_parameters[myChain,2], 
			alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, progressGraphs = FALSE, start_values = FALSE, zStart = zStart)
	}

	for(myChain in 1:nChains){
		kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
	}
	file_names <- list.files(outputDirs[1])	# the name of files are the same for each folder

	#	connections for saving the MCMC output corresponding to the target posterior distribution
	zConnection_target <- file(paste0(getwd(),"/zValues.txt"),open = "w")
	yConnection_target <- file(paste0(getwd(),"/yValues.txt"),open = "w")
	sigmainvConnection_target <- file(paste0(getwd(),"/sigmainvValues.txt"),open = "w")
	omegainvConnection_target <- file(paste0(getwd(),"/omegainvValues.txt"),open = "w")
	muConnection_target <- file(paste0(getwd(),"/muValues.txt"),open = "w")
	wConnection_target <- file(paste0(getwd(),"/wValues.txt"),open = "w")
	LambdaConnection_target <- vector('list',length = Kmax)
	cllConnection_target <- file(paste0(getwd(),"/cllValues.txt"),open = "w")
	for(k in 1:Kmax){
		LambdaConnection_target[[k]] <- file(paste0(getwd(),"/LambdaValues",k,".txt"),open = "w")  #K x p x q per iteration
	}




	
	#par(mfrow = c(1, 2))
	for( iteration in 2:mCycles ){
		
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			truncatedDirichletProcess(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 3, thinning = 1, burn = 2, alpha_prior= alpha_prior, g = g_h_parameters[myChain,1], h = g_h_parameters[myChain,2], 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, progressGraphs = FALSE, start_values = TRUE)
			kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
		}


		chains <- sample(nChains - 1, 1)
		chains <- c(chains, chains + 1)
		omegaINV1 <- as.numeric(read.table( paste0(outputDirs[ chains[1] ],'/omegainvValues.txt') ))
		omegaINV2 <- as.numeric(read.table( paste0(outputDirs[ chains[2] ],'/omegainvValues.txt') ))
		mh_denom <- sum( dgamma( omegaINV1, shape = g_h_parameters[chains[1],1], rate = g_h_parameters[chains[1],2], log = TRUE ) ) + 
				sum( dgamma( omegaINV2, shape = g_h_parameters[chains[2],1], rate = g_h_parameters[chains[2],2], log = TRUE ) )
		mh_nom   <- sum( dgamma( omegaINV1, shape = g_h_parameters[chains[2],1], rate = g_h_parameters[chains[2],2], log = TRUE )) + 
				sum( dgamma( omegaINV2, shape = g_h_parameters[chains[1],1], rate = g_h_parameters[chains[1],2], log = TRUE ) )
		mh_ratio <- mh_nom - mh_denom
		if( log(runif(1)) < mh_ratio ){
#		if( 0 < 1 ){

			# dir1 to tmp
			file.copy( 
					from      = paste0(outputDirs[ chains[1] ], '/', file_names), 
					to        = 'tmpDir',
					overwrite = TRUE
				)
			# dir2 to dir1
			file.copy( 
					from      = paste0(outputDirs[ chains[2] ], '/', file_names), 
					to        = outputDirs[ chains[1] ],
					overwrite = TRUE
				)
			# tmp to dir2
			file.copy( 
					from      = paste0('tmpDir', '/', file_names), 
					to        = outputDirs[ chains[2] ],
					overwrite = TRUE
				)
			
			mh_acceptance_rate <- mh_acceptance_rate + 1
		}

		for(myChain in 1:nChains){
			kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
		}

		z <- as.numeric(read.table('alpha_1/zValues.txt'))
		if( (iteration %% 50) == 0 ){
			par(mfrow = c(1,3))
			matplot(kValues[1:iteration, ], type = "l")
			points(1:iteration, kValues[1:iteration, 1], type = "b", col = 1)
			matplot(t(x_data), type = "l", col = as.numeric(as.factor(z)))
			matplot(t(originalX), type = "l", col = as.numeric(as.factor(z)))
			ar <- round(100*mh_acceptance_rate/iteration, 3)
			cat(paste0('-    iteration: ',iteration,'. Metropolis-Hastings acceptance rate: ', ar), '\n')
		}
		if(iteration %% thinning == 0){
			if(iteration > burnCycles){
				y        <- as.numeric(read.table('alpha_1/yValues.txt'))
				w        <- as.numeric(read.table('alpha_1/wValues.txt')) 
				mu       <- as.numeric(read.table('alpha_1/muValues.txt'))
				omegainv <- as.numeric(read.table('alpha_1/omegainvValues.txt'))
				sigmainv <- as.numeric(read.table('alpha_1/sigmainvValues.txt')) 
				cll      <- as.numeric(read.table('alpha_1/k.and.logl.Values.txt') )[2]
				for(k in 1:Kmax){
					Lambda <- as.numeric( read.table(paste0('alpha_1/LambdaValues', k, '.txt') ) )
					cat(Lambda, file = LambdaConnection_target[[k]], '\n', append = TRUE)
				}
				cat(z       , file = zConnection_target, '\n', append = TRUE)
				cat(y       , file = yConnection_target, '\n', append = TRUE)
				cat(w       , file = wConnection_target, '\n', append = TRUE)
				cat(mu      , file = muConnection_target, '\n', append = TRUE)
				cat(omegainv, file = omegainvConnection_target, '\n', append = TRUE)
				cat(sigmainv, file = sigmainvConnection_target, '\n', append = TRUE)
				cat(cll     , file = cllConnection_target, '\n', append = TRUE)
			}
		}
	}
	stopImplicitCluster()
	close(zConnection_target)
	close(yConnection_target)
	close(wConnection_target)
	close(muConnection_target)
	close(sigmainvConnection_target)
	close(omegainvConnection_target)
	close(cllConnection_target)
	for(k in 1:Kmax){
		close(LambdaConnection_target[[k]])
	}
	keepedSeq <- seq(burnCycles + thinning, mCycles, by = thinning)
	if( burnCycles > 0){
		burnedSeq <- 1:burnCycles
		write.table(file = 'burn_in_period_k.txt', kValues[burnedSeq, ])
	}
	write.table(file = 'kValues.txt', kValues[keepedSeq, ], quote = FALSE, row.names = FALSE, col.names = 1:nChains)
	setwd("../")
	cat('-    DONE.','\n')
}





