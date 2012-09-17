\name{OC.Isop}
\alias{OC.Isop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Iso P linie
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OC.Isop(p, n, k, resolution = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{resolution}{
%%     ~~Describe \code{resolution} here~~
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
function (p, n, k, resolution = 100) 
{
    sigma.max <- sigmap0(p)
    resolution <- 100
    sigmas <- seq(1e-06, sigma.max * (1 - 1e-05), length.out = resolution)
    mus <- rep(NA, resolution)
    integralwert <- mus
    for (i in 3:resolution) {
        mus[i] <- uniroot(mup, interval = c(0, 1e+06), tol = 1e-15, 
            sigma = sigmas[i], plevel = p)$root
    }
    for (i in 3:resolution) {
        integralwert[i] <- OC(mus[i], sigmas[i], n, k)
    }
    return(integralwert)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
