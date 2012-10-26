\name{quickOCIsop}
\alias{quickOCIsop}

\title{
Calculate a quick proxy for the complete OCIsop
}
\description{
Instead of using the complete iso-p-line this function evlauates only
for the maximum \code{sigma} possible for that p. More often than not
this yields the same aprroximative plan than using the whole line.
}
\usage{
quickOCIsop(p, n, k, type = "LR")
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
  \item{type}{
type must be either \"LR\" or \"BSK\"
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
\code{OCIsop}
}
\examples{
}

\keyword{ ~htest }

