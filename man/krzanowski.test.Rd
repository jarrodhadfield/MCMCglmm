\name{krzanowski.test}
\alias{krzanowski.test}
\title{Krzanowski's Comparison of Subspaces}
\description{Calculates statistics of Krzanowski's comparison of subspaces.}
\usage{
krzanowski.test(CA, CB, vecsA, vecsB, corr = FALSE, ...)
}
\arguments{
  \item{CA}{Matrix A}
  \item{CB}{Matrix B}
  \item{vecsA}{Vector of integers indexing the eigenvectors determining the 
subspace of A}
  \item{vecsB}{Vector of integers indexing the eigenvectors determining the 
subspace of B}
  \item{corr}{logical; if \code{TRUE} the variances of A and B are standardised}
  \item{...}{further arguments to be passed}
}
\value{
  \item{sumofS}{metric for overall similarity with 0 indicting no similarity and 
a value of \code{length(vecsA)} for identical subspaces}
  \item{angles}{angle in degrees between each best matched pair of vectors}
  \item{bisector}{vector that lies between each best matched pair of vectors}
}
\references{Krzanowski, W.J. (2000) Principles of Multivariate Analysis. OUP}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
CA<-rIW(diag(5),10, n=1)
CB<-rIW(diag(5),10, n=1)
krzanowski.test(CA, CB, vecsA=1:2, vecsB=1:2)
krzanowski.test(CA, CA, vecsA=1:2, vecsB=1:2)
}
\keyword{multivariate}

