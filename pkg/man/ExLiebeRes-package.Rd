\name{ExLiebeRes-package}
\alias{ExLiebeRes}
\docType{package}
\title{
Caculating optimal variable sampling plans of type \"LR\" or \"BSK\".
}
\description{
For both types of plans the OC is implemented. Using a four step
approximation generating intermediate plans the final optimal plan are
determined quickly. Details can be found in the reference resp. an
upcoming paper by the package author and Prof. Dr. Krumbholz. 
}
\details{
\tabular{ll}{
Package: \tab ExLiebeRes\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-10-22\cr
License: \tab GPL >= 2\cr
}

For an user the function \code{findOptPlan} is the main interface. For
any usable set of parameters p1, p2, alpha, beta the optimal variable
samplnig plan either
of type Lieberman/Reshnikov (\"LR\") or Bruhn-Suhr/Krumbholz (\"BSK\") 
is determinded.

For researchers the underlying OC function and intermediate plans may
also be evaluated directly.

Not all possible combinations of parameters give numerically usable
results. The underlying integrals become hard to evaluate if the sample
size gets big (n >= ca 750) or the difference between p1 and p2 very
small. Normally the internal \code{integrate} function just bails out in
these cases.


}
\author{
  The paper accompanies a paper written by Wolf Krumbholz and Detlef
  Steuer (to appear)
  
Maintainer: steuer@hsu-hh.de
}
\references{
paper to appear
}

\keyword{ package }
\seealso{

}
\examples{
findOptPlan(0.01, 0.06, 0.1, 0.1)
}
