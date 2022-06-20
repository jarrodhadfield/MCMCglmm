\name{knorm}
\alias{knorm}
\title{(Mixed) Central Moments of a Multivariate Normal Distribution}
\description{
  Forms a tensor of (mixed) central moments of a multivariate normal distribution
}
\usage{
knorm(V, k)
}
\arguments{
  \item{V}{(co)variance matrix}
  \item{k}{kth central moment, must be even}
}
\value{
   tensor
}
\references{Schott, J.R.(2003) Journal of Multivariate Analysis 87 (1) 177-190}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{dnorm}}}
\examples{
V<-diag(2)
knorm(V,2)
knorm(V,4)
}
\keyword{distribution}
