\name{kunif}
\alias{kunif}
\title{Central Moments of a Uniform Distribution}
\description{
 Returns the central moments of a uniform distribution
}
\usage{
kunif(min, max, k)
}
\arguments{
  \item{min, max}{lower and upper limits of the distribution.  Must be finite.}
  \item{k}{k central moment, must be even}
}
\value{
  kth central moment
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{dunif}}}
\examples{
kunif(-1,1,4)
y<-runif(1000,-1,1)
mean((y-mean(y))^4)
}
\keyword{distribution}

