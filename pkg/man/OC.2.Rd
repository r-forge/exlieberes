\name{OC.2}
\alias{OC.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ OC 2. Teil
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OC.2(mu, sigma, n, k, resolution = 2501)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (mu, sigma, n, k, resolution = 2501) 
{
    betap <- (n - 2)/2
    xgitter <- seq(0, qbeta(k, betap, betap), length.out = resolution)[-resolution]
    resolution <- resolution - 1
    werte <- rep(NA, resolution)
    gridwidth <- (xgitter[resolution] - xgitter[1])/resolution
    for (i in 1:resolution) {
        werte[i] <- innerOC.2(xgitter[i], mu, sigma, n, k)
    }
    sum(werte) * gridwidth
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
