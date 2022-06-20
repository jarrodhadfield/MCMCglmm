\name{commutation}
\alias{commutation}
\title{Commutation Matrix}
\description{
 Forms an mn x mn commutation matrix which transforms \eqn{vec({\bf A})}{vec(A)} into \eqn{vec({\bf A}^{'})}{vec(A')}, where \eqn{{\bf A}}{A} is an m x n matrix
}
\usage{
commutation(m, n)
}
\arguments{
  \item{m}{integer; number of rows of A}
  \item{n}{integer; number of columns of A}
}
\value{
  Commutation Matrix
}
\references{Magnus, J. R. & Neudecker, H. (1979) Annals of Statistics 7 (2) 381-394}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
commutation(2,2)
}
\keyword{array}

