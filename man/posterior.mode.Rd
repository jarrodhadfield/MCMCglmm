\name{posterior.mode}
\alias{posterior.mode}
\title{Estimates the marginal parameter modes using kernel density estimation}
\description{Estimates the marginal parameter modes using kernel density estimation}
\usage{
posterior.mode(x, adjust=0.1, ...)
}
\arguments{
 \item{x}{mcmc object}
 \item{adjust}{numeric, passed to \code{\link{density}} to adjust the bandwidth of the kernal density}
 \item{...}{other arguments to be passed}
}
\value{
  modes of the kernel density estimates
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{density}}}
\examples{
v<-rIW(as.matrix(1),10, n=1000)
hist(v)
abline(v=posterior.mode(mcmc(v)), col="red")
}
\keyword{distribution}

