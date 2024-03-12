\name{sir}
\alias{sir}
\title{Design Matrix for Simultaneous and Recursive Relationships between Responses}
\description{Forms design matrix for simultaneous and recursive relationships between responses}
\usage{
sir(formula1=NULL, formula2=NULL, diag0=FALSE)
}
\arguments{
  \item{formula1}{formula}
  \item{formula2}{formula}
  \item{diag0}{logical: should the design matrix have zero's along the diagonal}
}
\value{
    design matrix 
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
fac1<-factor(sample(letters[1:3], 5, TRUE), levels=letters[1:3])
fac2<-factor(sample(letters[1:3], 5, TRUE), levels=letters[1:3])
cbind(fac1, fac2)
sir(~fac1, ~fac2)

}
\keyword{manip}

