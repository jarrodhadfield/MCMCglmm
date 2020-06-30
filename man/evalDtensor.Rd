\name{evalDtensor}
\alias{evalDtensor}
\title{Evaluates a list of (mixed) partial derivatives}
\description{
 Evaluates a list of (mixed) partial derivatives
}
\usage{
evalDtensor(x, mu, m=1)
}
\arguments{
  \item{x}{unevaluated (list) of expression(s)}
  \item{mu}{values at which the derivatives are evaluated: names need to match 
terms in x}
  \item{m}{order of derivative}
}
\value{
   tensor 
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{Dtensor}}, \code{\link{D}}}
\examples{
f<-expression(beta_1 + time*beta_2+u)
Df<-Dtensor(f, eval=FALSE, m=2)
evalDtensor(Df, mu=data.frame(beta_1=0.5, beta_2=1, time=3, u=2.3))
Dtensor(f, mu=c(1,3,1,2.3), m=2)
}
\keyword{math}
\keyword{array}
