\name{mult.memb}
\alias{mult.memb}
\title{Design Matrices for Multiple Membership Models}
\description{Forms design matrices for multiple membership models}
\usage{
mult.memb(formula)
}
\arguments{
  \item{formula}{formula}
}
\value{
    design matrix
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\details{
Currently \code{mult.memb} can only usefully be used inside an \code{idv} variance function. The formula usually contains serveral factors that have the same factor levels.
}
\examples{
fac1<-factor(sample(letters[1:3], 5, TRUE), levels=letters[1:3])
fac2<-factor(sample(letters[1:3], 5, TRUE), levels=letters[1:3])
cbind(fac1, fac2)
mult.memb(~fac1+fac2)
}
\keyword{manip}

