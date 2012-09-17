\name{integrandOC.1}
\alias{integrandOC.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ erster Integrand OC
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
integrandOC.1(t, n, deltaL, deltaU, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{deltaL}{
%%     ~~Describe \code{deltaL} here~~
}
  \item{deltaU}{
%%     ~~Describe \code{deltaU} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (t, n, deltaL, deltaU, x) 
{
    wt <- sqrt(t * (n - 1))
    dchisq(t, n - 1) * (pnorm(-deltaU + (2 * x - 1) * wt) - pnorm(-deltaL + 
        wt))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
