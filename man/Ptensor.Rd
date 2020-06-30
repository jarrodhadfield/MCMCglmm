\name{Ptensor}
\alias{Ptensor}
\title{Tensor of Sample (Mixed) Central Moments}
\description{Forms a tensor of sample (mixed) central moments}
\usage{
Ptensor(x, k)
}
\arguments{
  \item{x}{matrix; traits in columns samples in rows}
  \item{k}{kth central moment}
}
\value{
   \item{}{tensor}
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
n<-1000
y<-matrix(rnorm(n), n/2, 2)
Ptensor(y,2)
cov(y)*((n-1)/n)
}
\keyword{array}

