\name{posterior.ante}
\alias{posterior.ante}
\title{Posterior distribution of ante-dependence parameters}
\description{Posterior distribution of ante-dependence parameters}
\usage{
posterior.ante(x, vtype="ante1")
}
\arguments{
 \item{x}{mcmc object of (co)variances stacked column-wise}
 \item{vtype}{type of ante variance function used}
}
\value{
  posterior ante-dependence parameters (innovation variances followed by regression ceofficients)
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{posterior.cor}}, \code{\link{posterior.evals}}, \code{\link{posterior.inverse}}}
\examples{
v<-rIW(diag(2),10, n=1000)
plot(posterior.ante(mcmc(v),vtype="ante1"))
}
\keyword{distribution}

