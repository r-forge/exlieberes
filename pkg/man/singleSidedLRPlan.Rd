\name{singleSidedLRPlan}
\alias{singleSidedLRPlan}
\title{
Calculating an exact Lieberman/Reshnikov plan for a one-sided
specification limit.
}
\description{
Acoording to the given set of parameters p1, p2, alpha, beta the LR plan
for a single sided specification limit is calculated.
}
\usage{
singleSidedLRPlan(p1, p2, alpha, beta, verbosity.level = 0)
}
\arguments{
  \item{p1}{ requested AQL
}
  \item{p2}{
requested RQL
}
  \item{alpha}{
error of first kind
}
  \item{beta}{
error of second kind
}
  \item{verbosity.level}{
level verboseness: 0 no output only a plan structure returned, 1
generated plan and intermediate approximative plans get printed, >=2 detailed information about the search process
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{  A plan structure is returned as a list containing these items:
  \item{ntilde}{single sided LR plan sample size $n$ , }
  \item{ltilde}{critical value untransformed , }
  \item{ktilde}{critical value transformed}
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
singleSidedLRPlan(0.01, 0.06, 0.1, 0.1, verbosity.level=0)
singleSidedLRPlan(0.01, 0.06, 0.1, 0.1, verbosity.level=2)
## same result, but more information about intermediate steps
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ ~htest }

