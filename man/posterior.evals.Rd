\name{posterior.evals}
\alias{posterior.evals}
\title{Posterior distribution of eigenvalues}
\description{Posterior distribution of eigenvalues}
\usage{
posterior.evals(x)
}
\arguments{
 \item{x}{mcmc object of (co)variances stacked column-wise}
}
\value{
  posterior eigenvalues
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{posterior.cor}}, \code{\link{posterior.inverse}}, \code{\link{posterior.ante}}}
\examples{
v<-rIW(diag(2),3, n=1000)
hist(posterior.evals(mcmc(v))[,2])
}
\keyword{distribution}

