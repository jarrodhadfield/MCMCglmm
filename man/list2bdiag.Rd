\name{list2bdiag}
\alias{list2bdiag}
\title{Forms the direct sum from a list of matrices}
\description{Forms a block-diagonal matrix from a list of matrices}
\usage{
list2bdiag(x)
}
\arguments{
  \item{x}{list of square matrices}
}
\value{
    matrix
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
M<-list(rIW(diag(3), 10), rIW(diag(2), 10))
list2bdiag(M)
}
\keyword{manip}

