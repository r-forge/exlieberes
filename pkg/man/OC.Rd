\name{OC}
\alias{OC}
\title{
Calculate the operation characteric (OC) of a plan
}
\description{
For both types of plans \"LR\" and \"BSK\" the OC for a given (fitting)
pair of ( \code{mu}, \code{sigma}) is evaluated.
}
\usage{
OC(mu, sigma, n, k, resolution = 100, type = "LR")
}
\arguments{
  \item{mu}{
   \code{mu}, is a mu(sigma), i.e. must match with sigma, see reference.
}
  \item{sigma}{
 \code{sigma}
}
  \item{n}{ samplesize
}
  \item{k}{ critical value
}
  \item{resolution}{
defines the resolution for the internal numerical integrations
}
  \item{type}{
must be \"LR\" oder \"BSK\".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ real number, value of OC at that point.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Detlef Steuer, steuer@hsu-hh.de, http://fawn.hsu-hh.de/~steuer
}
\note{
Keep in mind that not all combinations of parameters give
results. Integration may easily fail, if \code{mu} and \code{sigma} do
not fit. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ ~htest }

