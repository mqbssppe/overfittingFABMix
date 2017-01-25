\name{update_z4}
\alias{update_z4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update_z4(w, mu, Lambda, SigmaINV, K, x_data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{Lambda}{
%%     ~~Describe \code{Lambda} here~~
}
  \item{SigmaINV}{
%%     ~~Describe \code{SigmaINV} here~~
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
function (w, mu, Lambda, SigmaINV, K, x_data) 
{
    n <- dim(x_data)[1]
    p <- dim(x_data)[2]
    probs <- array(data = 0, dim = c(n, K))
    for (k in 1:K) {
        center_x <- x_data - matrix(mu[k, ], nrow = n, ncol = p, 
            byrow = TRUE)
        x_var <- Lambda[k, , ] \%*\% t(Lambda[k, , ])
        diag(x_var) <- diag(x_var) + 1/diag(SigmaINV)
        lpdf <- log(w[k]) + dmvnorm(center_x, mean = rep(0, p), 
            sigma = x_var, log = TRUE)
        probs[, k] <- lpdf
    }
    probs <- array(t(apply(probs, 1, function(tmp) {
        return(exp(tmp - max(tmp)))
    })), dim = c(n, K))
    z <- apply(probs, 1, function(tmp) {
        return(sample(K, 1, prob = tmp))
    })
    results <- vector("list", length = 2)
    names(results) <- c("w", "z")
    results[[1]] <- w
    results[[2]] <- z
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
