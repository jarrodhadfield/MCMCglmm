\name{Dtensor}
\alias{Dtensor}
\title{Tensor of (mixed) partial derivatives}
\description{Forms tensor of (mixed) partial derivatives}
\usage{
Dtensor(expr, name=NULL, mu = NULL, m=1, evaluate = TRUE)
}
\arguments{
  \item{expr}{'expression'}
  \item{name}{character vector, giving the variable names with respect to which 
derivatives will be computed. If NULL all variables in the expression will be 
used}
  \item{mu}{optional: numeric vector, at which the derivatives are evaluated}
  \item{m}{order of derivative}
  \item{evaluate}{logical; if \code{TRUE} the derivatives are evaluated at 
\code{mu}, if \code{FALSE} the derivatives are left unevaluated}
}
\value{
  \item{Dtensor}{(list) of unevaluated expression(s) if \code{evaluate=FALSE} or 
a tensor if \code{evaluate=TRUE}}
}
\references{Rice, S.H. (2004) Evolutionary Theory: Mathematical and Conceptual 
Foundations. Sinauer (MA) USA.}
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\seealso{\code{\link{evalDtensor}}, \code{\link{Dexpressions}}, \code{\link{D}}}
\examples{
f<-expression(beta_1 + time * beta_2 + u)
Dtensor(f,eval=FALSE)
}
\keyword{math}
\keyword{array}

