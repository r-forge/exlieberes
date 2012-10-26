\name{findOptPlan}
\alias{findOptPlan}
\title{
Determine the optimal variables sampling plan of given type
}
\description{
For a given type of plan the optimal variable sampling plan is determined.
}
\usage{
findOptPlan(p1, p2, alpha, beta, type = "LR", use.quickApproxPlan = TRUE, given.plan = NULL, verbosity.level = 0, show.plots = FALSE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{type}{
type must be either \"LR\" or \"BSK\"
}
  \item{use.quickApproxPlan}{
if TRUE a \code{quickApproxPlan} is used to shorten search time. Always
faster. FALSE only possible for research.
}
  \item{given.plan}{
a plan structure may be inserted in \code{findOptPlan} to evaluate it or
to use is a start plan for further optimisation.
}
  \item{verbosity.level}{
level verboseness: 0 no output only a plan structure returned, 1
generetad plan gets printed, 2 detailed information about the search process
}
  \item{show.plots}{
if TRUE the OC is plotted for intermediate steps.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ a plan structure }
\references{
upcoming paper
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~htest}

