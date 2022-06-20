\name{pkk}
\alias{pkk}
\title{Probability that all multinomial categories have a non-zero count.}
\description{
  Calculates the probability that all categories in a multinomial have a non-zero count.
  }
\usage{
pkk(prob, size)
}
\arguments{
  \item{prob}{numeric non-negative vector of length K, specifying the probability for the K classes; is internally normalized to sum 1. Infinite and missing values are not allowed.}
  \item{size}{integer, say N, specifying the total number of objects that are put into K boxes in the typical multinomial experiment.}
}
\value{
  probability that there is at least one object in each of the K boxes
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
p<-runif(4)
pkk(p, 10)
}
\keyword{manip}


