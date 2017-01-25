\name{compute_sufficient_statistics}
\alias{compute_sufficient_statistics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_sufficient_statistics(y, z, K, x_data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{x_data}{
%%     ~~Describe \code{x_data} here~~
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
function (y, z, K, x_data) 
{
    cluster_size <- numeric(K)
    p <- dim(x_data)[2]
    q <- dim(y)[2]
    sx <- array(data = 0, dim = c(K, p))
    sy <- array(data = 0, dim = c(K, q))
    sxx <- array(data = 0, dim = c(K, p, p))
    syy <- array(data = 0, dim = c(K, q, q))
    sxy <- array(data = 0, dim = c(K, p, q))
    for (k in 1:K) {
        index <- which(z == k)
        cluster_size[k] <- length(index)
        if (cluster_size[k] > 0) {
            sx[k, ] <- colSums(array(x_data[index, ], dim = c(cluster_size[k], 
                p)))
            sy[k, ] <- colSums(array(y[index, ], dim = c(cluster_size[k], 
                q)))
            for (i in index) {
                sxx[k, , ] <- sxx[k, , ] + x_data[i, ] \%*\% t(x_data[i, 
                  ])
                syy[k, , ] <- syy[k, , ] + y[i, ] \%*\% t(y[i, 
                  ])
                sxy[k, , ] <- sxy[k, , ] + x_data[i, ] \%*\% t(y[i, 
                  ])
            }
        }
    }
    results <- vector("list", length = 6)
    names(results) <- c("cluster_size", "sx", "sy", "sxx", "syy", 
        "sxy")
    results[[1]] <- cluster_size
    results[[2]] <- sx
    results[[3]] <- sy
    results[[4]] <- sxx
    results[[5]] <- syy
    results[[6]] <- sxy
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
