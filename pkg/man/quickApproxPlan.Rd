\name{quickApproxPlan}
\alias{quickApproxPlan}
\title{
quickApproxPlan
}
\description{
an approximation for the looked after plan is generated without
calculating the complete OC function. Instead only one point of the OC for the
maximum value of \code{sigma} is evalutated and used as a surprisingly not
so crude approximation. Obviously this gains great benefits in computing
time. The result is used as an intermediate before searchin in detail
using the complete OC
}
\usage{
quickApproxPlan(p1, p2, alpha, beta, type = "LR", verbosity.level = 0)
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
  \item{type}{
is either "LR" (Leibermann/Reshnikov) or "BSK" (Bruhn-suhr/Krumbholz)
indicating the type of plan that shall be generated.
}
  \item{verbosity.level}{
level verboseness: 0 no output only a plan structure returned, 1
generetad plan gets printed, 2 detailed information about the search process
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
  Detlef Steuer, steuer@hsu-hh.de, http://fawn.hsu-hh.de/~steuer
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
function (p1, p2, alpha, beta, type = "LR", verbosity.level = 0) 
{
    if (!((type == "LR") | (type == "BSK"))) 
        stop("Type must be either \"LR\" or \"BSK\"!")
    startplan <- singleSidedLRPlan(p1, p2, alpha, beta, verbosity.level = verbosity.level)
    nmin <- startplan$ntilde
    if (type == "LR") {
        kopt <- startplan$ktilde
    }
    if (type == "BSK") {
        kopt <- pnorm(startplan$ltilde/sqrt(startplan$ntilde))
    }
    kdiff <- 0.02 * kopt
    PlanGefunden <- FALSE
    direction <- "unknown"
    while (!PlanGefunden) {
        if (verbosity.level >= 2) 
            cat("Testing plan p1=", p1, " p2 = ", p2, " n = ", 
                nmin, " k =", kopt, "\n")
        Lminp1 <- quickOCIsop(p1, nmin, kopt, type = type)
        if (verbosity.level >= 2) 
            cat("L_min(p1) =", Lminp1)
        Lmaxp2 <- quickOCIsop(p2, nmin, kopt, type = type)
        if (verbosity.level >= 2) 
            cat(" L_max(p2) =", Lmaxp2, "\n")
        if ((Lminp1 >= 1 - alpha) & (Lmaxp2 <= beta)) {
            PlanGefunden <- TRUE
            plan <- list(ntilde = nmin, ltilde = NA, ktilde = kopt)
        }
        else if (!(Lminp1 >= 1 - alpha) & !(Lmaxp2 < beta)) {
            nmin <- nmin + 1
            direction <- "unknown"
            if (verbosity.level >= 2) 
                cat("increase n now: ", nmin, "\n")
        }
        else {
            if (Lminp1 < 1 - alpha) {
                if (direction == "unknown") 
                  direction <- "up"
                if (direction == "up") 
                  kopt <- kopt + kdiff
                else {
                  direction <- "up"
                  kdiff <- kdiff/2
                  kopt <- kopt + kdiff
                }
                if (verbosity.level >= 2) 
                  cat("new k_opt: ", kopt, "\n")
            }
            else {
                if (direction == "unknown") 
                  direction <- "down"
                if (direction == "down") 
                  kopt <- kopt - kdiff
                else {
                  direction <- "down"
                  kdiff <- kdiff/2
                  kopt <- kopt - kdiff
                }
                if (verbosity.level >= 2) 
                  cat("new k_opt: ", kopt, "\n")
            }
        }
    }
    if (verbosity.level >= 1) 
        cat("Fast approximative plan for p_1 = ", p1, " p_2 = ", 
            p2, " alpha = ", alpha, " beta = ", beta, ": n = ", 
            plan$ntilde, " k = ", plan$ktilde, "\n")
    return(plan)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
