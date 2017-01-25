\name{update_all_y}
\alias{update_all_y}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update_all_y(x_data, mu, SigmaINV, Lambda, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_data}{
%%     ~~Describe \code{x_data} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{SigmaINV}{
%%     ~~Describe \code{SigmaINV} here~~
}
  \item{Lambda}{
%%     ~~Describe \code{Lambda} here~~
}
  \item{z}{
%%     ~~Describe \code{z} here~~
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
function (x_data, mu, SigmaINV, Lambda, z) 
{
    p <- dim(Lambda)[2]
    q <- dim(Lambda)[3]
    n <- dim(x_data)[1]
    n_k <- table(z)
    alive <- as.numeric(names(n_k))
    y_new <- array(data = 0, dim = c(n, q))
    for (k in alive) {
        ind <- which(z == k)
        center_x <- x_data[ind, ] - matrix(mu[k, ], nrow = as.numeric(n_k[as.character(k)]), 
            ncol = p, byrow = TRUE)
        Alpha <- t(Lambda[k, , ]) \%*\% SigmaINV \%*\% Lambda[k, 
            , ]
        diag(Alpha) <- diag(Alpha) + 1
        Alpha <- solve(Alpha)
        tmp <- Alpha \%*\% t(Lambda[k, , ]) \%*\% SigmaINV
        y_mean <- t(apply(center_x, 1, function(tk) {
            return(tmp \%*\% tk)
        }))
        if (q == 1) {
            y_mean <- t(y_mean)
        }
        y_new[ind, ] <- t(apply(y_mean, 1, function(tk) {
            return(mvrnorm(n = 1, mu = tk, Sigma = Alpha))
        }))
    }
    results <- vector("list", length = 1)
    results[[1]] <- y_new
    names(results) <- c("y")
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
