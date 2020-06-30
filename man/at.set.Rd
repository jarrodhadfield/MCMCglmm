\name{at.set}
\alias{at.set}
\title{Incidence Matrix of Combined Levels within a Factor}
\description{Incidence Matrix of Combined Levels within a Factor}
\usage{
at.set(x, level)
}
\arguments{
  \item{x}{factor}
  \item{level}{set of factor levels}
}
\value{
    incidence matrix for the set \code{level} in x
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{at.level}}}
\examples{
fac<-gl(3,10,30, labels=letters[1:3])
x<-rnorm(30)
model.matrix(~at.set(fac,2:3):x)
}
\keyword{manip}

