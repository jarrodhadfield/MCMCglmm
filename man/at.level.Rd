\name{at.level}
\alias{at.level}
\title{Incidence Matrix of Levels within a Factor}
\description{Incidence matrix of levels within a factor}
\usage{
at.level(x, level)
}
\arguments{
  \item{x}{factor}
  \item{level}{factor level}
}
\value{
    incidence matrix for level in x
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{at.set}}}
\examples{
fac<-gl(3,10,30, labels=letters[1:3])
x<-rnorm(30)
model.matrix(~at.level(fac,"b"):x)
}
\keyword{manip}

