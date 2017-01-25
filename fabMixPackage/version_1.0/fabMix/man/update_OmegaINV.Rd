\name{update_OmegaINV}
\alias{update_OmegaINV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update_OmegaINV(Lambda, K, g, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lambda}{
%%     ~~Describe \code{Lambda} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
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
function (Lambda, K, g, h) 
{
    p <- dim(Lambda)[2]
    q <- dim(Lambda)[3]
    OmegaINV <- array(data = 0, dim = c(q, q))
    betaVector <- numeric(q)
    for (k in 1:K) {
        betaVector <- betaVector + colSums(array(Lambda[k, , 
            ]^2, dim = c(p, q)))
    }
    diag(OmegaINV) <- rgamma(q, shape = g + K * p/2, rate = h + 
        betaVector/2)
    return(OmegaINV)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
