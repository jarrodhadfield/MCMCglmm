\name{Ddivergence}
\alias{Ddivergence}
\title{d-divergence}
\description{Calculates Ovaskainen's (2008) d-divergence between 2 zero-mean multivariate normal distributions.}
\usage{
Ddivergence(CA=NULL, CB=NULL, n=10000)
}
\arguments{
  \item{CA}{Matrix A}
  \item{CB}{Matrix B}
  \item{n}{number of Monte Carlo samples for approximating the integral}
}
\value{
  d-divergence
}
\note{In versions of MCMCglmm <2.26 Ovaskainen's (2008) d-divergence was incorrectly calculated.
}	
\references{Ovaskainen, O.  et. al. (2008) Proc. Roy. Soc - B (275) 1635 593-750}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
CA<-rIW(diag(2),10, n=1)
CB<-rIW(diag(2),10, n=1)
Ddivergence(CA, CB)
}
\keyword{multivariate}
