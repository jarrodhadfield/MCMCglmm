\name{posterior.ante}
\alias{posterior.ante}
\title{Posterior distribution of ante-dependence parameters}
\description{Posterior distribution of ante-dependence parameters}
\usage{
posterior.ante(x,k=1)
}
\arguments{
 \item{x}{mcmc object of (co)variances stacked column-wise}
 \item{k}{order of the ante-dependence structure}
}
\value{
  posterior ante-dependence parameters (innovation variances followed by regression ceofficients)
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{posterior.cor}}, \code{\link{posterior.evals}}, \code{\link{posterior.inverse}}}
\examples{
v<-rIW(diag(2),10, n=1000)
plot(posterior.ante(mcmc(v),1))
}
\keyword{distribution}

