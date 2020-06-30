\name{posterior.cor}
\alias{posterior.cor}
\title{Transforms posterior distribution of covariances into correlations}
\description{Transforms posterior distribution of covariances into correlations}
\usage{
posterior.cor(x)
}
\arguments{
 \item{x}{mcmc object of (co)variances stacked column-wise}
}
\value{
  posterior correlation matrices
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{posterior.evals}}, \code{\link{posterior.inverse}}, \code{\link{posterior.ante}}}
\examples{
v<-rIW(diag(2),3, n=1000)
hist(posterior.cor(mcmc(v))[,2])
}
\keyword{distribution}

