\name{OCIsop}
\alias{OCIsop}
\title{
Calculate OC along an iso-p-line.
}
\description{
Along an iso-p-line for the given p the OC is evaluated on an
equidistant grid.
}
\usage{
OCIsop(p, n, k, resolution = 100, type = "LR")
}
\arguments{
  \item{p}{
probability for which the iso-p-line is defined
}
  \item{n}{
sample size
}
  \item{k}{
critical value
}
  \item{resolution}{
number of equidistant points on the iso-p-line-
}
  \item{type}{
type must be either \"LR\" or \"BSK\"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
array of flaoting point numbers of length resolution
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~htest }

