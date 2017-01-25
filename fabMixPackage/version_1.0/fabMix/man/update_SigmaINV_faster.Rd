\name{update_SigmaINV_faster}
\alias{update_SigmaINV_faster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update_SigmaINV_faster(x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_data}{
%%     ~~Describe \code{x_data} here~~
}
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{Lambda}{
%%     ~~Describe \code{Lambda} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{alpha_sigma}{
%%     ~~Describe \code{alpha_sigma} here~~
}
  \item{beta_sigma}{
%%     ~~Describe \code{beta_sigma} here~~
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
function (x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma) 
{
    p <- dim(Lambda)[2]
    q <- dim(Lambda)[3]
    n <- length(z)
    SigmaINV <- array(data = 0, dim = c(p, p))
    s <- numeric(p)
    alive <- as.numeric(names(table(z)))
    for (k in alive) {
        ind <- which(z == k)
        n_k <- length(ind)
        tmp <- matrix(mu[k, ], nrow = n_k, ncol = p, byrow = TRUE) + 
            t(apply(array(y[ind, ], dim = c(n_k, q)), 1, function(tk) {
                return(array(Lambda[k, , ], dim = c(p, q)) \%*\% 
                  tk)
            }))
        s <- s + colSums((x_data[ind, ] - tmp)^2)
    }
    diag(SigmaINV) <- rgamma(p, shape = alpha_sigma + n/2, rate = beta_sigma + 
        s/2)
    return(SigmaINV)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
