\name{Tri2M}
\alias{Tri2M}
\title{Lower/Upper Triangle Elements of a Matrix}
\description{
  Lower/Upper triangle elements of a matrix or forms a matrix from a vector of 
lower/upper triangle elements
  }
\usage{
Tri2M(x, lower.tri = TRUE, reverse = TRUE, diag = TRUE)
}
\arguments{
  \item{x}{Matrix or vector}
  \item{lower.tri}{If \code{x} is a matrix then the lower triangle (\code{TRUE}) 
or upper triangle \code{FALSE} elements (including diagonal elements) are 
returned. If \code{x} is a vector a matrix is formed under the assumption that 
\code{x} are the lower triangle (\code{TRUE}) or upper triangle (\code{FALSE}) 
elements.}
  \item{reverse}{logical: if\code{TRUE} a symmetric matrix is formed, if 
\code{FALSE} the remaining triangle is left as zeros.}
  \item{diag}{logical: if\code{TRUE} diagonal elements are included.}
}
\value{
  \item{}{numeric or matrix}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
M<-rIW(diag(3), 10)
x<-Tri2M(M)
x
Tri2M(x, reverse=TRUE)
Tri2M(x, reverse=FALSE)
}
\keyword{manip}


