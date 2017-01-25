\name{dealWithLabelSwitching_same_sigma}
\alias{dealWithLabelSwitching_same_sigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dealWithLabelSwitching_same_sigma(x_data, outputFolder, q, burn, z.true, compute_regularized_expression, Km)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_data}{
%%     ~~Describe \code{x_data} here~~
}
  \item{outputFolder}{
%%     ~~Describe \code{outputFolder} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{burn}{
%%     ~~Describe \code{burn} here~~
}
  \item{z.true}{
%%     ~~Describe \code{z.true} here~~
}
  \item{compute_regularized_expression}{
%%     ~~Describe \code{compute_regularized_expression} here~~
}
  \item{Km}{
%%     ~~Describe \code{Km} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x_data, outputFolder, q, burn, z.true, compute_regularized_expression, 
    Km) 
{
    p <- dim(x_data)[2]
    n <- dim(x_data)[1]
    if (missing(Km)) {
        Km <- 20
    }
    if (missing(burn)) {
        burn <- 1:1
    }
    if (missing(compute_regularized_expression)) {
        compute_regularized_expression = FALSE
    }
    cat(paste0("-    (5) Dealing with label switching for q = ", 
        q), "\n")
    setwd(outputFolder)
    cat(paste0("         * Entering directory: ", getwd()), "\n")
    z <- as.matrix(read.table("zValues.txt"))
    logl <- read.table("kValues.txt", header = T)
    logl <- logl[-burn, ]
    z <- z[-burn, ]
    K <- as.numeric(names(sort(table(logl[, 1]), decreasing = TRUE)[1]))
    kSelected <- K
    cat(paste0("         * Posterior mode corresponds to K = ", 
        K), "\n")
    if (K == 1) {
        cat(paste0("         *  no label switching algorithms are applied"), 
            "\n")
    }
    else {
        index <- which(logl[, 1] == K)
        Kindex <- index
        logl <- logl[index, 2]
        z <- z[index, ]
        mapIndex <- which(logl == max(logl))[1]
        zPivot <- as.numeric(z[mapIndex, ])
        K <- max(z)
        skiniko <- FALSE
        if (K < Km) {
            skiniko <- TRUE
            kLAST <- K
            newZ <- t(apply(z, 1, function(y) {
                myI <- which(y == K)
                y[myI] <- rep(Km, length(myI))
                return(y)
            }))
            z <- newZ
            newZ <- zPivot
            myI <- which(newZ == K)
            zPivot[myI] <- rep(Km, length(myI))
            K <- max(z)
        }
        m <- length(logl)
        if (missing(z.true)) {
            ls <- label.switching(method = c("ECR", "ECR-ITERATIVE-1"), 
                zpivot = zPivot, z = z, K = K)
        }
        else {
            ls <- label.switching(method = c("ECR", "ECR-ITERATIVE-1"), 
                zpivot = zPivot, z = z, K = K, groundTruth = z.true)
        }
        oldLabels <- 1:K
        index <- Kindex
        allocationsECR <- z
        for (i in 1:m) {
            myPerm <- order(ls$permutations$ECR[i, ])
            allocationsECR[i, ] <- myPerm[z[i, ]]
        }
        write.table(allocationsECR, file = "reordered_allocations_ecr.txt")
        cat(paste0("         * write file: `reordered_allocations_ecr.txt`"), 
            "\n")
        l <- as.matrix(read.table(paste0("LambdaValues", 1, ".txt")))
        J <- dim(l)[2]
        mcmc <- array(data = NA, dim = c(m, K, J))
        for (k in 1:K) {
            l <- as.matrix(read.table(paste0("LambdaValues", 
                k, ".txt")))
            l <- l[-burn, ]
            mcmc[, k, ] <- l[index, ]
        }
        if (skiniko == TRUE) {
            tmp1 <- mcmc[, kLAST, ]
            tmp2 <- mcmc[, Km, ]
            mcmc[, kLAST, ] <- tmp2
            mcmc[, Km, ] <- tmp1
        }
        lambda.perm.mcmc <- permute.mcmc(mcmc, ls$permutations$ECR)$output
        write.table(lambda.perm.mcmc, file = "reordered_lambda_ecr.txt")
        cat(paste0("         * write file: `reordered_lambda_ecr.txt`"), 
            "\n")
        lambda.mean <- array(data = NA, dim = c(K, p, q))
        lambda.map <- array(data = NA, dim = c(K, p, q))
        for (k in 1:K) {
            lMean <- apply(lambda.perm.mcmc[, k, ], 2, mean)
            lambda.mean[k, , ] <- matrix(lMean, nrow = p, ncol = q, 
                byrow = TRUE)
            lambda.map[k, , ] <- matrix(lambda.perm.mcmc[mapIndex, 
                k, ], nrow = p, ncol = q, byrow = TRUE)
        }
        write.table(lambda.mean, file = "lambda_estimate_ecr.txt", 
            col.names = paste("lambda", 1:p, rep(1:q, each = p), 
                sep = "_"))
        cat(paste0("         * write file: `lambda_estimate_ecr.txt`"), 
            "\n")
        mu <- read.table("muValues.txt")[max(burn) + Kindex, 
            ]
        mu.mcmc <- array(data = NA, dim = c(m, K, p))
        for (k in 1:K) {
            mu.mcmc[, k, ] <- as.matrix(mu[, k + K * ((1:p) - 
                1)])
        }
        if (skiniko == TRUE) {
            tmp1 <- mu.mcmc[, kLAST, ]
            tmp2 <- mu.mcmc[, Km, ]
            mu.mcmc[, kLAST, ] <- tmp2
            mu.mcmc[, Km, ] <- tmp1
        }
        mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations$ECR)$output
        write.table(mu.mcmc, file = "reordered_mu_ecr.txt")
        cat(paste0("         * write file: `reordered_mu_ecr.txt`"), 
            "\n")
        mu.mean <- array(data = NA, dim = c(K, p))
        mu.map <- array(data = NA, dim = c(K, p))
        for (k in 1:K) {
            for (j in 1:p) {
                mu.mean[k, j] <- mean(mu.mcmc[, k, j])
                mu.map[k, j] <- mu.mcmc[mapIndex, k, j]
            }
        }
        write.table(mu.mean, file = "mu_estimate_ecr.txt")
        cat(paste0("         * write file: `mu_estimate_ecr.txt`"), 
            "\n")
        w.mcmc <- array(as.matrix(read.table("wValues.txt")[max(burn) + 
            Kindex, ]), dim = c(length(Kindex), K, 1))
        w.mcmc_raw <- w.mcmc
        w.mcmc <- array(w.mcmc[, 1:K, ], dim = c(length(Kindex), 
            K, 1))
        if (skiniko == TRUE) {
            tmp1 <- w.mcmc[, kLAST, ]
            tmp2 <- w.mcmc[, Km, ]
            w.mcmc[, kLAST, ] <- tmp2
            w.mcmc[, Km, ] <- tmp1
        }
        w.mcmc <- permute.mcmc(w.mcmc, ls$permutations$ECR)$output
        w.mean <- numeric(K)
        w.map <- numeric(K)
        for (k in 1:K) {
            w.mean[k] <- mean(w.mcmc[, k, 1])
            w.map[k] <- w.mcmc[mapIndex, k, 1]
        }
        write.table(w.mcmc, file = "reordered_weights_ecr.txt")
        cat(paste0("         * write file: `reordered_weights_ecr.txt`"), 
            "\n")
        write.table(w.mean, file = "weights_estimate_ecr.txt")
        cat(paste0("         * write file: `weights_estimate_ecr.txt`"), 
            "\n")
        write.table(ls$clusters, file = "singleBestClusterings.txt", 
            quote = FALSE, row.names = FALSE)
        cat(paste0("         * write file: `singleBestClusterings.txt`"), 
            "\n")
        zMAP <- as.matrix(read.table("singleBestClusterings.txt", 
            header = TRUE))[1, ]
        mapAllocationsPosteriorProbs <- numeric(n)
        for (i in 1:n) {
            mapAllocationsPosteriorProbs[i] <- length(which(allocationsECR[, 
                i] == zMAP[i]))
        }
        mapAllocationsPosteriorProbs <- mapAllocationsPosteriorProbs/dim(allocationsECR)[1]
        write.table(mapAllocationsPosteriorProbs, file = "classificationProbabilities.txt")
        cat(paste0("         * write file: `classificationProbabilities.txt`"), 
            "\n")
        aliveClusters <- as.numeric(names(table(ls$clusters[1, 
            ])))
        sigmaINV <- as.matrix(read.table("sigmainvValues.txt")[Kindex, 
            ])
        covmat <- array(data = 0, dim = c(length(aliveClusters), 
            p, p))
        rownames(covmat) <- as.character(aliveClusters)
        for (iter in 1:m) {
            for (k in aliveClusters) {
                lambda_k <- matrix(lambda.perm.mcmc[iter, k, 
                  ], nrow = p, ncol = q, byrow = T)
                covmat[as.character(k), , ] <- covmat[as.character(k), 
                  , ] + lambda_k \%*\% t(lambda_k)
                diag(covmat[as.character(k), , ]) <- diag(covmat[as.character(k), 
                  , ]) + as.numeric(1/sigmaINV[iter, ])
            }
        }
        covmat <- covmat/m
        for (k in 1:length(aliveClusters)) {
            write.table(covmat[k, , ], file = paste0("estimated_cov_cluster_", 
                k, ".txt"))
            cat(paste0("         * write file: `estimated_cov_cluster_", 
                k, ".txt`"), "\n")
            write.table(lambda.map[k, , ], file = paste0("lambda_map_", 
                k, ".txt"))
            cat(paste0("         * write file: `lambda_map_", 
                k, ".txt`"), "\n")
        }
        if (compute_regularized_expression == TRUE) {
            cat(paste("-    computing regularized expressions..."), 
                "\n")
            yValues <- read.table("yValues.txt")[Kindex, ]
            regularizedExpression <- array(data = 0, dim = c(length(aliveClusters), 
                p))
            regularizedExpression2 <- array(data = 0, dim = c(length(aliveClusters), 
                p, q))
            rownames(regularizedExpression) <- rownames(regularizedExpression2) <- as.character(aliveClusters)
            for (iter in 1:m) {
                for (k in aliveClusters) {
                  lambda_k <- matrix(lambda.perm.mcmc[iter, k, 
                    ], nrow = p, ncol = q, byrow = T)
                  yMean <- array(data = 0, dim = c(q, 1))
                  index_k <- as.numeric(which(allocationsECR[iter, 
                    ] == k))
                  for (i in index_k) {
                    yMean <- yMean + t(yValues[iter, (1:q - 1) * 
                      n + i])
                  }
                  yMean <- yMean/length(index_k)
                  tmp <- t(apply(lambda_k, 1, function(y) {
                    y * yMean
                  }))
                  regularizedExpression2[as.character(k), , ] <- regularizedExpression2[as.character(k), 
                    , ] + tmp
                  regularizedExpression[as.character(k), ] <- regularizedExpression[as.character(k), 
                    ] + rowSums(tmp)
                }
                if (iter\%\%50 == 0) {
                  cat(paste("iter = ", iter), "\n")
                }
            }
            regularizedExpression <- regularizedExpression/m
            regularizedExpression2 <- regularizedExpression2/m
            for (j in 1:q) {
                write.table(regularizedExpression2[, , j], file = paste0("estimated_regularized_expression_per_cluster_", 
                  j, ".txt"))
                cat(paste0("         * write file: `estimated_regularized_expression_per_cluster_", 
                  j, ".txt`"), "\n")
            }
            write.table(regularizedExpression, file = "estimated_regularized_expression_per_cluster.txt")
            cat(paste0("         * write file: `estimated_regularized_expression_per_cluster.txt`"), 
                "\n")
        }
        cat(paste0("-    Done."), "\n")
    }
    setwd("../")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
