\name{sm2asreml}
\alias{sm2asreml}
\title{Converts sparseMatrix to asreml's giv format}
\description{
 Converts sparseMatrix to asreml's giv format: row-ordered, upper triangle 
sparse matrix. 
}
\usage{
sm2asreml(A=NULL, rownames=NULL)
}
\arguments{
  \item{A}{sparseMatrix}
  \item{rownames}{rownames of A}
}
\value{
  data.frame: if \code{A} was formed from a pedigree equivalent to giv format 
returned by \code{asreml.Ainverse}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{inverseA}
\examples{
data(bird.families)
A<-inverseA(bird.families)
Aasreml<-sm2asreml(A$Ainv, A$node.names)
}
\keyword{manip}


