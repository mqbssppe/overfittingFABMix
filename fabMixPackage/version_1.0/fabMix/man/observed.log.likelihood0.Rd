\name{observed.log.likelihood0}
\alias{observed.log.likelihood0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
observed.log.likelihood0(x_data, w, mu, Lambda, Sigma, z)
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
  \item{Sigma}{
%%     ~~Describe \code{Sigma} here~~
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
function (x_data, w, mu, Lambda, Sigma, z) 
{
    p <- dim(x_data)[2]
    n <- dim(x_data)[1]
    ct <- -(p/2) * log(2 * pi)
    probs <- numeric(n)
    if (is.null(z)) {
        alive <- 1
    }
    else {
        alive <- as.numeric(names(table(z)))
    }
    newW <- numeric(length(w))
    newW[alive] <- w[alive]/sum(w[alive])
    loggedValues <- array(data = NA, dim = c(n, length(alive)))
    colnames(loggedValues) <- alive
    for (k in alive) {
        center_x <- x_data - matrix(mu[k, ], nrow = n, ncol = p, 
            byrow = TRUE)
        x_var <- Lambda[k, , ] \%*\% t(Lambda[k, , ])
        diag(x_var) <- diag(x_var) + Sigma
        loggedValues[, as.character(k)] <- log(newW[k]) + dmvnorm(center_x, 
            mean = rep(0, p), sigma = x_var, log = TRUE)
    }
    lMax <- apply(loggedValues, 1, max)
    if (length(alive) == 1) {
        logL <- sum(lMax + log(exp(apply(loggedValues, 2, function(y) {
            return(y - lMax)
        }))))
    }
    else {
        logL <- sum(lMax + log(rowSums(exp(apply(loggedValues, 
            2, function(y) {
                return(y - lMax)
            })))))
    }
    return(logL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
