\name{sigmap0}
\alias{sigmap0}
\title{sigmap0()
}
\description{
For a given $p < 0.5$ the function returns an upper limit on $sigma$
for calculating ISO-p lines.
}
\usage{
sigmap0(p)
}
\arguments{
  \item{p}{
p - numerical value 0 < p < 0.5 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a float, the upper limit on sigma for calculating
the OC along an ISO-p line for the given $p$.
}
\references{
Krumbholz/Bruhn-Suhr 1990, Statistical papers 31 (p. 195-207), equation 1.18
}
\author{
Detlef Steuer, steuer@hsu-hh.de
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
sigmap0(0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

