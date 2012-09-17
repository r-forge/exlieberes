\name{innerOC.2}
\alias{innerOC.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ innere OC
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
innerOC.2(x, mu, sigma, n, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
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
function (x, mu, sigma, n, k) 
{
    dL <- deltaL(mu, sigma, n)
    if (DEBUG) 
        cat("x ", x, "mu ", mu, "sigma ", sigma, "n ", n, "k ", 
            k, "\n")
    result <- integrate(ininOC.2, 0, grossA(psi(x, k, (n - 2)/2), 
        x, sigma, n), subdivisions = 2500, x = x, deltaL = dL, 
        n = n)$value
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
