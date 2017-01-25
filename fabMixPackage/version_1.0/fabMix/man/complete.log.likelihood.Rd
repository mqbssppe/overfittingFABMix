\name{complete.log.likelihood}
\alias{complete.log.likelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
complete.log.likelihood(x_data, w, mu, Lambda, SigmaINV, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_data}{
%%     ~~Describe \code{x_data} here~~
}
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
function (x_data, w, mu, Lambda, SigmaINV, z) 
{
    n <- length(z)
    p <- dim(Lambda)[2]
    q <- dim(Lambda)[3]
    probs <- numeric(n)
    alive <- as.numeric(names(table(z)))
    for (k in alive) {
        index <- which(z == k)
        center_x <- x_data[index, ] - matrix(mu[k, ], nrow = length(index), 
            ncol = p, byrow = TRUE)
        x_var <- Lambda[k, , ] \%*\% t(Lambda[k, , ])
        diag(x_var) <- diag(x_var) + 1/diag(SigmaINV)
        x_var <- try(solve(x_var), TRUE)
        if (is.numeric(x_var) == TRUE) {
            probs[index] <- log(w[k]) - 0.5 * apply(center_x, 
                1, function(tmp) {
                  return(as.numeric(t(tmp) \%*\% x_var \%*\% tmp))
                }) + 0.5 * log(det(x_var))
        }
    }
    return(sum(probs))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
