\name{fOptPlan}
\alias{fOptPlan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Alte langsame Funktion
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fOptPlan(p1, p2, alpha, beta)
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
    startplan <- LRP(p1, p2, alpha, beta)
    nmin <- startplan$ntilde
    kopt <- startplan$ktilde
    kdiff <- 2e-05
    PlanGefunden <- FALSE
    direction <- "unknown"
    while (!PlanGefunden) {
        cat("Teste Plan mit p1=", p1, " p2= ", p2, " n= ", nmin, 
            " und k=", kopt, "\n")
        result1 <- OC.Isop(p1, nmin, kopt)
        plot(result1[result1 > 0.1])
        Lminp1 <- min(result1[result1 > 0.1], na.rm = TRUE)
        cat("Lminp1 =", Lminp1)
        result2 <- OC.Isop(p2, nmin, kopt)
        plot(result2[result2 > 0.02])
        Lmaxp2 <- max(result2, na.rm = TRUE)
        cat(" Lmaxp2 =", Lmaxp2, "\n")
        if ((Lminp1 >= 1 - alpha) & (Lmaxp2 <= beta)) {
            PlanGefunden <- TRUE
        }
        else if (!(Lminp1 >= 1 - alpha) & !(Lmaxp2 < beta)) {
            nmin <- nmin + 1
            direction <- "unknown"
            cat("n\303\244chstes n: ", nmin, "\n")
        }
        else {
            if (Lminp1 < 1 - alpha) {
                if (direction == "up") 
                  kopt <- kopt + kdiff
                else {
                  kdiff <- kdiff/2
                  kopt <- kopt + kdiff
                  direction <- "up"
                }
                cat("neues kopt: ", kopt, "\n")
            }
            else {
                if (direction == "down") 
                  kopt <- kopt - kdiff
                else {
                  kdiff <- kdiff/2
                  kopt <- kopt - kdiff
                  direction <- "down"
                }
                cat("neues kopt: ", kopt, "\n")
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
