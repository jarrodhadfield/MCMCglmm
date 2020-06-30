\name{Dexpressions}
\alias{Dexpressions}
\title{List of unevaluated expressions for (mixed) partial derivatives of 
fitness with respect to linear predictors.}
\description{
  Unevaluated expressions for (mixed) partial derivatives of fitness with 
respect to linear predictors for survival and fecundity.
}
\usage{
  Dexpressions
}
\value{
  \item{PW.d0W}{Fitness (W) function for the Poisson-Weibull (PW) model.}
  \item{PW.d1Wds}{First Partial derivative of fitness (d1W) with respect to 
survival (d1s) linear predictor for the Poisson-Weibull (PW) model.}
  \item{PW.d1Wdf}{First Partial derivative of fitness (d1W) with respect to 
fecundity (d1f) linear redictor for the Poisson-Weibull (PW) model.} 
  \item{PW.d3Wd2sd1f}{Mixed third partial derivative of fitness (d3W) with 2nd 
derivative of survival linear predictor (d2s) and first derivative of fecundity 
linear predictor (d1f) from the Poisson-Weibull (PW) model.}   
  \item{PW.d3Wdsd2f}{and so on ...}
  \item{PW.d2Wd2f}{}
  \item{PW.d2Wd2s}{}
  \item{PW.d3Wd3s}{}
  \item{PW.d3Wd3f}{}
  }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{Dtensor}}}
\keyword{symbolmath}


