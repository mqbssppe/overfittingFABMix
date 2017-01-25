\name{getStuffForDIC}
\alias{getStuffForDIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getStuffForDIC(x_data, outputFolder, q, burn, Km, normalize, discardLower)
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
  \item{Km}{
%%     ~~Describe \code{Km} here~~
}
  \item{normalize}{
%%     ~~Describe \code{normalize} here~~
}
  \item{discardLower}{
%%     ~~Describe \code{discardLower} here~~
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
function (x_data, outputFolder, q, burn, Km, normalize, discardLower) 
{
    cat(paste0("-    (4) Computing information criteria for q = ", 
        q), "\n")
    if (missing(normalize)) {
        normalize = TRUE
    }
    if (normalize) {
        x_data <- scale(x_data, center = TRUE, scale = TRUE)
        cat("-    NOTE: using standardizing data.", "\n")
    }
    n <- dim(x_data)[1]
    p <- dim(x_data)[2]
    if (missing(Km)) {
        Km <- 20
    }
    if (missing(burn)) {
        burn <- 1:1
    }
    setwd(outputFolder)
    cat(paste0("         - Entering directory: ", getwd()), "\n")
    z <- as.matrix(read.table("zValues.txt"))
    logl <- read.table("cllValues.txt")
    tmp <- apply(z, 1, function(y) {
        length(table(y))
    })
    logl <- cbind(tmp, logl)
    z <- z[-burn, ]
    logl <- logl[-burn, ]
    cat(paste0("            Nclusters:    ", paste(as.character(names(table(logl[, 
        1]))), collapse = "     ")), "\n")
    cat(paste0("            Frequency:    ", paste(as.character(as.numeric(table(logl[, 
        1]))), collapse = "    ")), "\n")
    K <- as.numeric(names(sort(table(logl[, 1]), decreasing = TRUE)[1]))
    kSelected <- K
    index <- which(logl[, 1] == K)
    Kindex <- index
    m <- length(index)
    ECR <- matrix(1:Km, nrow = m, ncol = Km, byrow = T)
    permutations <- vector("list", length = 1)
    permutations[[1]] <- ECR
    names(permutations) <- "ECR"
    ls <- vector("list", length = 1)
    names(ls) <- "permutations"
    ls$permutations <- permutations
    l <- as.matrix(read.table(paste0("LambdaValues", 1, ".txt")))
    J <- dim(l)[2]
    mcmc <- array(data = NA, dim = c(m, Km, J))
    for (k in 1:Km) {
        l <- as.matrix(read.table(paste0("LambdaValues", k, ".txt")))
        l <- l[-burn, ]
        mcmc[, k, ] <- l[index, ]
    }
    lambda.perm.mcmc <- permute.mcmc(mcmc, ls$permutations$ECR)$output
    for (k in 1:Km) {
        lMean <- apply(lambda.perm.mcmc[, k, ], 2, mean)
    }
    mu <- read.table("muValues.txt")[max(burn) + Kindex, ]
    mu.mcmc <- array(data = NA, dim = c(m, Km, p))
    for (k in 1:Km) {
        mu.mcmc[, k, ] <- as.matrix(mu[, k + Km * ((1:p) - 1)])
    }
    mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations$ECR)$output
    SigmaINV <- as.matrix(read.table("sigmainvValues.txt")[max(burn) + 
        Kindex, ])
    SigmaINV.mcmc <- SigmaINV
    Sigma.mcmc <- 1/SigmaINV.mcmc
    w.mcmc <- array(as.matrix(read.table("wValues.txt")[max(burn) + 
        Kindex, ]), dim = c(length(Kindex), Km, 1))
    w.mcmc <- array(w.mcmc[, 1:Km, ], dim = c(length(Kindex), 
        Km, 1))
    w.mcmc <- permute.mcmc(w.mcmc, ls$permutations$ECR)$output
    lValues <- numeric(m)
    cll <- 0
    aic <- 0
    bic <- 0
    maxL <- 1
    i <- 1
    lambda.current <- array(data = NA, dim = c(Km, p, q))
    Sigma.current <- Sigma.mcmc[i, ]
    mu.current <- mu.mcmc[i, , ]
    if (Km == 1) {
        Sigma.current <- array(Sigma.current, dim = c(1, p))
        mu.current <- array(mu.mcmc[i, , ], dim = c(1, p))
    }
    for (k in 1:Km) {
        ldraw <- lambda.perm.mcmc[i, k, ]
        lambda.current[k, , ] <- matrix(ldraw, nrow = p, ncol = q, 
            byrow = TRUE)
        for (i1 in 1:p) {
            if (Sigma.current[i1] > 100) {
                cat(paste0("oops: ", i), "\n")
                Sigma.current[i1] <- 100
            }
        }
    }
    obsL <- observed.log.likelihood0(x_data = x_data, w = w.mcmc[i, 
        , 1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, 
        z = z[i, ])
    lValues[i] <- obsL
    maxL <- obsL
    cll <- cll + obsL
    aic <- aic + obsL
    bic <- bic + obsL
    iterMax <- i
    for (i in 2:m) {
        lambda.current <- array(data = NA, dim = c(Km, p, q))
        Sigma.current <- Sigma.mcmc[i, ]
        mu.current <- mu.mcmc[i, , ]
        if (Km == 1) {
            Sigma.current <- array(Sigma.current, dim = c(1, 
                p))
            mu.current <- array(mu.mcmc[i, , ], dim = c(1, p))
        }
        for (k in 1:Km) {
            ldraw <- lambda.perm.mcmc[i, k, ]
            lambda.current[k, , ] <- matrix(ldraw, nrow = p, 
                ncol = q, byrow = TRUE)
            for (i1 in 1:p) {
                if (Sigma.current[i1] > 100) {
                  cat(paste0("oops: ", i), "\n")
                  Sigma.current[i1] <- 100
                }
            }
        }
        obsL <- observed.log.likelihood0(x_data = x_data, w = w.mcmc[i, 
            , 1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, 
            z = z[i, ])
        lValues[i] <- obsL
        if (obsL > maxL) {
            maxL <- obsL
            iterMax <- i
        }
        cll <- cll + obsL
        aic <- aic + obsL
        bic <- bic + obsL
    }
    if (missing(discardLower)) {
        discardLower <- 0.1
    }
    if (discardLower == FALSE) {
        cll <- cll/m
    }
    else {
        cll <- mean(lValues[which(lValues > as.numeric(quantile(lValues, 
            discardLower)))])
    }
    dic_classicMAP <- -4 * cll + 2 * maxL
    dic_starMAP <- -6 * cll + 4 * maxL
    dic_classic <- dic_classicMAP
    dic_star <- dic_starMAP
    aic <- -2 * aic/m + 2 * (kSelected * (p + p * q - q * (q - 
        1)/2) + p + kSelected - 1)
    bic <- -2 * bic/m + log(n) * (kSelected * (p + p * q - q * 
        (q - 1)/2) + p + kSelected - 1)
    aic_MAX <- -2 * maxL + 2 * (kSelected * (p + p * q - q * 
        (q - 1)/2) + p + kSelected - 1)
    bic_MAX <- -2 * maxL + log(n) * (kSelected * (p + p * q - 
        q * (q - 1)/2) + p + kSelected - 1)
    dic <- c(aic, bic, dic_classic, dic_star, dic_classicMAP, 
        dic_starMAP, aic_MAX, bic_MAX)
    names(dic) <- c("AIC", "BIC", "DIC1", "DIC*2", "DIC", "DIC_2", 
        "AIC_map", "BIC_map")
    write.table(file = "informationCriteria_map_model.txt", dic[c(5, 
        6, 7, 8)], col.names = paste0("q_", q), quote = FALSE)
    write.table(file = "lValues_map.txt", lValues, quote = FALSE)
    setwd("../")
    cat(paste0("         - Information criteria written to `", 
        outputFolder, "/informationCriteria_map_model.txt`."), 
        "\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
