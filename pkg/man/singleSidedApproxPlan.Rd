\name{singleSidedApproxPlan}
\alias{singleSidedApproxPlan}
\title{
  Calculate a first approximating plan for a one-sided LR plan
}
\description{
An approximative plan for the LR plan based on the approximation of Uhlmann
generated. This plan may generate an
$n$ minimally larger than the one of the exact one-sided LR plan.
}
\usage{
singleSidedApproxPlan(p1, p2, alpha, beta, verbosity.level = 0)
}
\arguments{
  \item{p1}{
AQL
}
  \item{p2}{
RQL
}
  \item{alpha}{
error of first kind
}
  \item{beta}{
error of second kind
}
  \item{verbosity.level}{
level of verboseness: 0 no output besides plan structure returned, >= 1
some information about the search process printed
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ A plan structure is returned as a list containing these items:
  \item{ntilde}{aprroximation of sample size $n$}
  \item{ltilde}{approximation of critical value untransformed}
  \item{ktilde}{approximation of critical value transformed}
}
\references{
Uhlmann, Statistische Qualitätskontrolle, 2. Aufl., Stuttgart (1982),
p. 186 ff
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
singleSidedApproxPlan(0.01, 0.06, 0.1, 0.1)
}
\keyword{ htest }

