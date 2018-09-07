
obs.log.likelihood0 <- function(x_data, w, mu, Lambda, Sigma, z){
        p <- dim(x_data)[2]
        n <- dim(x_data)[1]
        ct <- -(p/2)*log(2*pi)
	K <- length(w)
        probs <- numeric(n)
	alive <- 1:K
        newW <- numeric(length(w))
        newW[alive] <- w[alive]/sum(w[alive])
        loggedValues <- array(data = NA, dim = c(n, length(alive)))
        colnames(loggedValues) <- alive
        for(k in alive){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
                x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
                diag(x_var) <- diag(x_var) + Sigma
                loggedValues[ ,as.character(k)] <- log(newW[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
        }
        lMax <- apply(loggedValues, 1, max)
        if( length(alive) == 1 ){
                logL <- sum(lMax + log( exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ) ) )
        }else{
                logL <- sum(lMax + log( rowSums(exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ))))
        }
        return( logL )
}


compute_standard_BIC <- function(sameSigma = TRUE, x_data, outputFolder, q, burn, Km, normalize, discardLower){
        cat(paste0('-    (4) Computing information criteria for q = ', q), '\n')
        if(missing(normalize)){normalize = TRUE}
        if(normalize){
                x_data <- scale(x_data, center = TRUE, scale = TRUE)
                cat('-    NOTE: using standardizing data.','\n')
        }
        n <- dim(x_data)[1]
        p <- dim(x_data)[2]
        if(q == 0){sameSigma = FALSE}
        if(missing(Km)){Km <- 0}
        if(missing(burn)){burn <- 0}
        setwd(outputFolder)
        cat(paste0('         - Entering directory: ', getwd()),'\n')
        z <- as.matrix(read.table("zValues.txt"))
        tmp  <- apply(z,1,function(y){length(table(y))})
        if(burn > 0){
                z <- z[-(1:burn),]
        }
        K <- Km
        kSelected <- K
	m <- dim(z)[1]
        index <- 1:m
        Kindex <- index

        #this is artificial

        ECR <- matrix(1:Km, nrow = m, ncol = Km, byrow=T)
        permutations <- vector('list', length = 1)
        permutations[[1]] <- ECR
        names(permutations) <- "ECR"
        ls <- vector('list', length = 1)
        names(ls) <- "permutations"
        ls$permutations <- permutations 
        #Lambda
        if(q > 0){
                l <- as.matrix(read.table(paste0("LambdaValues",1,".txt")))
                J <- dim(l)[2]
                mcmc <- array(data = NA, dim = c(m,Km,J))
                for(k in 1:Km){
                        l <- as.matrix(read.table(paste0("LambdaValues",k,".txt")))
                        if(burn > 0){
                        l <- l[-(1:burn),]}
                        mcmc[,k,] <- l[index,]
                }
                lambda.perm.mcmc <- permute.mcmc(mcmc, ls$permutations$ECR)$output
                for(k in 1:Km){
                        lMean <- apply(lambda.perm.mcmc[,k,],2,mean)
                }
        }
        mu <- read.table("muValues.txt")# auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
        if(burn > 0){
                mu <- mu[-(1:burn),] 
        }
        mu <- mu[Kindex,]
        mu.mcmc <- array(data = NA, dim = c(m,Km,p))
        for(k in 1:Km){
                mu.mcmc[,k,] <- as.matrix(mu[,k + Km*((1:p)-1)])
        }
        mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations$ECR)$output
        #
        if(sameSigma == TRUE){
                SigmaINV <- as.matrix(read.table("sigmainvValues.txt"))# auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
                if(burn > 0){
                        SigmaINV <- SigmaINV[-(1:burn),] 
                }
                SigmaINV <- SigmaINV[Kindex, ] 
                SigmaINV.mcmc <- SigmaINV
        }else{
                SigmaINV <- as.matrix(read.table("sigmainvValues.txt")) # auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
                if(burn > 0){
                        SigmaINV <- SigmaINV[-(1:burn),] 
                }
                SigmaINV <- SigmaINV[Kindex,  ] 
                SigmaINV.mcmc <- array(data = NA, dim = c(m,Km,p))
                for(k in 1:Km){
                        SigmaINV.mcmc[,k,] <- as.matrix(SigmaINV[,((k-1)*p + 1):(k*p)])
                }
                SigmaINV.mcmc <- permute.mcmc(SigmaINV.mcmc, ls$permutations$ECR)$output

        }
        Sigma.mcmc <- 1/SigmaINV.mcmc
        #SigmaINV.mean <- as.numeric(apply(SigmaINV,2,mean))
        w.mcmc <- as.matrix(read.table("wValues.txt"))
        w.mcmc <- array(w.mcmc, dim = c(dim(w.mcmc)[1], Km, 1))
        if(burn > 0){
                w.mcmc <- w.mcmc[-(1:burn),,]
                w.mcmc <- w.mcmc[Kindex,]
        }else{
                w.mcmc <- w.mcmc[Kindex,,]
        }
        w.mcmc <- array(w.mcmc[,1:Km],dim = c(length(Kindex),Km,1))
        w.mcmc <- permute.mcmc(w.mcmc, ls$permutations$"ECR")$output

        lValues <- numeric(m)
        cll <- 0
        aic <- 0
        bic <- 0
        maxL <- 1
        i <- 1
        if(sameSigma == TRUE){
                Sigma.current <- Sigma.mcmc[i, ]
        }else{
                Sigma.current <- Sigma.mcmc[i, , ]              
        }
        mu.current <- mu.mcmc[i,,]
        if( Km == 1 ){  
                Sigma.current <- array(Sigma.current, dim = c(1, p))
                mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
        }
        if(q > 0){lambda.current <- array(data = NA, dim = c(Km,p,q))}
        for(k in 1:Km){
                if(q > 0){
                        ldraw <- lambda.perm.mcmc[i,k,]
                        lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
                }
                for(i1 in 1:p){
                        if( Sigma.current[ i1] > 100 ){ cat(paste0('oops: ', i),'\n'); Sigma.current[i1] <- 100 }
                }
        }
        if(sameSigma == TRUE){
                obsL <- obs.log.likelihood0(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
        }else{
                if(q > 0){
                        obsL <- obs.log.likelihood0_Sj(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
                }else{
                        obsL <- obs.log.likelihood0_Sj_q0(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Sigma = Sigma.current, z = z[i,])
                }
        }
        lValues[i] <- obsL
        maxL <- obsL
        cll <- cll + obsL
        aic <- aic + obsL
        bic <- bic + obsL
        iterMax <- i
        if(sameSigma == TRUE){
                for(i in 2:m){
        #               cat(paste0("i  = ", i), "\n")
                        lambda.current <- array(data = NA, dim = c(Km,p,q))
                        Sigma.current <- Sigma.mcmc[i, ]
                        mu.current <- mu.mcmc[i,,]
                        if( Km == 1 ){  
                                Sigma.current <- array(Sigma.current, dim = c(1, p))
                                mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
                        }
                        for(k in 1:Km){
        #                       cat(paste0("  k  = ", k), "\n")
                                ldraw <- lambda.perm.mcmc[i,k,]
                                lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
                                for(i1 in 1:p){
                                        if( Sigma.current[ i1] > 100 ){ cat(paste0('oops: ', i),'\n'); Sigma.current[ i1] <- 100 }
                                }
                        }
                        obsL <- obs.log.likelihood0(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
                        lValues[i] <- obsL
                        if( obsL > maxL ){
                                maxL <- obsL    
                                iterMax <- i
                        }
                        cll <- cll + obsL
                        aic <- aic + obsL
                        bic <- bic + obsL
                }
        }else{
                for(i in 2:m){
        #               cat(paste0("i  = ", i), "\n")
                        lambda.current <- array(data = NA, dim = c(Km,p,q))
                        Sigma.current <- Sigma.mcmc[i, , ]
                        mu.current <- mu.mcmc[i,,]
                        if( Km == 1 ){  
                                Sigma.current <- array(Sigma.current, dim = c(1, p))
                                mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
                        }
                        for(k in 1:Km){
        #                       cat(paste0("  k  = ", k), "\n")
                                if(q > 0){
                                        ldraw <- lambda.perm.mcmc[i,k,]
                                        lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
                                }
                                for(i1 in 1:p){
                                        if( Sigma.current[k, i1] > 100 ){ Sigma.current[k, i1] <- 100 }
                                }
                        }
                        if(q > 0){
                                obsL <- obs.log.likelihood0_Sj(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
                        }else{
                                obsL <- obs.log.likelihood0_Sj_q0(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Sigma = Sigma.current, z = z[i,])
                        }
                        lValues[i] <- obsL
                        if( obsL > maxL ){
                                maxL <- obsL    
                                iterMax <- i
                        }
                        cll <- cll + obsL
                        aic <- aic + obsL
                        bic <- bic + obsL
                }
        }
        if(missing(discardLower)){ discardLower <- 0.1 }
        if ( discardLower == FALSE){
                cll <- cll/m
        }else{
                cll <- mean( lValues[which( lValues > as.numeric(quantile(lValues, discardLower)) )] )
        }
        dic_classicMAP <- -4*cll + 2*maxL
        dic_starMAP <- -6*cll + 4*maxL
        dic_classic <- dic_classicMAP
        dic_star <- dic_starMAP
        if(sameSigma == TRUE){
                aic <- -2*aic/m + 2*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 ) 
                bic <- -2*bic/m + log(n)*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 )
                aic_MAX <- -2*maxL + 2*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 ) 
                bic_MAX <- -2*maxL + log(n)*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 )
        }else{
                aic <- -2*aic/m + 2*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 ) 
                bic <- -2*bic/m + log(n)*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 )
                aic_MAX <- -2*maxL + 2*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 ) 
                bic_MAX <- -2*maxL + log(n)*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 )
        }
        dic <- c(aic, bic, dic_classic, dic_star, dic_classicMAP, dic_starMAP, aic_MAX, bic_MAX)
        names(dic) <- c('AIC', 'BIC', 'DIC1', 'DIC*2', 'DIC', 'DIC_2', 'AIC_map', 'BIC_map')    
#       write.table(file = 'informationCriteria_map_model.txt', dic[c(1,2,5,6,7,8)], col.names = paste0('q_',q), quote = FALSE)
        write.table(file = 'informationCriteria_map_model.txt', dic[c(5,6,7,8)], col.names = paste0('q_',q), quote = FALSE)
        write.table(file = 'lValues_map.txt', lValues, quote = FALSE)
        setwd("../")
        cat(paste0('         - Information criteria written to `', outputFolder,'/informationCriteria_map_model.txt`.'), '\n')
}


