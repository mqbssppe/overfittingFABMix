\name{compute_A_B_G_D_and_simulate_mu_Lambda}
\alias{compute_A_B_G_D_and_simulate_mu_Lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_A_B_G_D_and_simulate_mu_Lambda(SigmaINV, suff_statistics, OmegaINV, K, priorConst1, T_INV, v_r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SigmaINV}{
%%     ~~Describe \code{SigmaINV} here~~
}
  \item{suff_statistics}{
%%     ~~Describe \code{suff_statistics} here~~
}
  \item{OmegaINV}{
%%     ~~Describe \code{OmegaINV} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{priorConst1}{
%%     ~~Describe \code{priorConst1} here~~
}
  \item{T_INV}{
%%     ~~Describe \code{T_INV} here~~
}
  \item{v_r}{
%%     ~~Describe \code{v_r} here~~
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
function (SigmaINV, suff_statistics, OmegaINV, K, priorConst1, 
    T_INV, v_r) 
{
    p <- dim(suff_statistics$sx)[2]
    q <- dim(suff_statistics$sy)[2]
    A <- array(data = 0, dim = c(K, p, p))
    B <- mu <- array(data = 0, dim = c(K, p))
    G <- array(data = 0, dim = c(K, p, q, q))
    D <- array(data = 0, dim = c(K, p, q))
    mu <- array(data = 0, dim = c(K, p))
    Lambdas <- array(data = 0, dim = c(K, p, q))
    for (k in 1:K) {
        diag(A[k, , ]) <- 1/(suff_statistics$cluster_size[k] * 
            diag(SigmaINV) + diag(T_INV))
        for (r in 1:p) {
            subIndex <- 1:v_r[r]
            tmp <- OmegaINV[subIndex, subIndex] + suff_statistics$syy[k, 
                subIndex, subIndex] * SigmaINV[r, r] - A[k, r, 
                r] * (SigmaINV[r, r]^2) * suff_statistics$sy[k, 
                subIndex] \%*\% t(suff_statistics$sy[k, subIndex])
            G[k, r, subIndex, subIndex] <- solve(tmp)
            D[k, r, subIndex] <- suff_statistics$sxy[k, r, subIndex] * 
                SigmaINV[r, r] - A[k, r, r] * (SigmaINV[r, r]^2) * 
                suff_statistics$sx[k, r] * suff_statistics$sy[k, 
                subIndex]
            lambda_mean <- G[k, r, subIndex, subIndex] \%*\% D[k, 
                r, subIndex]
            Lambdas[k, r, subIndex] <- mvrnorm(n = 1, mu = lambda_mean, 
                Sigma = G[k, r, subIndex, subIndex])
        }
        B[k, ] <- SigmaINV \%*\% (suff_statistics$sx[k, ] - array(Lambdas[k, 
            , ], dim = c(p, q)) \%*\% array(suff_statistics$sy[k, 
            ], dim = c(q, 1))) + priorConst1
        mu_mean <- A[k, , ] \%*\% B[k, ]
        mu[k, ] <- mvrnorm(n = 1, mu = mu_mean, Sigma = A[k, 
            , ])
    }
    results <- vector("list", length = 6)
    results[[1]] <- A
    results[[2]] <- B
    results[[3]] <- G
    results[[4]] <- D
    results[[5]] <- Lambdas
    results[[6]] <- mu
    names(results) <- c("A", "B", "G", "D", "Lambdas", "mu")
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
