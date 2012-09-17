\name{LRP}
\alias{LRP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LRP
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LRP(p1, p2, alpha, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
%%     ~~Describe \code{p1} here~~
}
  \item{p2}{
%%     ~~Describe \code{p2} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
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
function (p1, p2, alpha, beta) 
{
    n1 <- ceiling(((qnorm(1 - alpha) - qnorm(beta))/(qnorm(p2) - 
        qnorm(p1)))^2)
    l1n <- qnorm(1 - alpha) + sqrt(n1) * qnorm(p1)
    l2n <- qnorm(beta) + sqrt(n1) * qnorm(p2)
    l1 <- (l1n + l2n)/2
    n0 <- floor(n1 * (1 + l1^2/(2 * n1)))
    l <- l1 * sqrt(1 + l1^2/(2 * n1))
    cat("N\303\244herungsplan n0=", n0, " l0=", l, "\n")
    return(list(ntilde = n0 - 1, ltilde = l, ktilde = pbeta(0.5 + 
        l/(2 * (n0 - 1)), (n0 - 2)/2, (n0 - 2)/2)))
    n <- n0 - 2
    while (n < n0 + 3) {
        cat("n= ", n, "\n")
        cat("(i)  ", pt(l, n - 1, sqrt(n) * qnorm(p1)), ">= 1-alpha", 
            1 - alpha, "\n")
        cat("(ii) ", pt(l, n - 1, sqrt(n) * qnorm(p2)), "<  beta", 
            beta, "\n")
        if (pt(l, n - 1, sqrt(n) * qnorm(p1)) >= (1 - alpha) & 
            pt(l, n - 1, sqrt(n) * qnorm(p2)) < beta) {
            cat("min n =", n, "\n")
            plan <- list(ntilde = n, ltilde = l, ktilde = pbeta(0.5 + 
                l/(2 * (n - 1)), (n - 2)/2, (n - 2)/2))
            return(plan)
        }
        n <- n + 1
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
