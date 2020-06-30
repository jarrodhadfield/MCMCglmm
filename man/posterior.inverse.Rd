\name{posterior.inverse}
\alias{posterior.inverse}
\title{Posterior distribution of matrix inverse}
\description{Posterior distribution of matrix inverse}
\usage{
posterior.inverse(x)
}
\arguments{
 \item{x}{mcmc object of (co)variances stacked column-wise}
}
\value{
  posterior of inverse (co)variance matrices
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{posterior.cor}}, \code{\link{posterior.evals}}, \code{\link{posterior.ante}}}
\examples{
v<-rIW(diag(2),3, n=1000)
plot(posterior.inverse(mcmc(v)))
}
\keyword{distribution}

