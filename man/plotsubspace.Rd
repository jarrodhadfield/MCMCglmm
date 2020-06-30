\name{plotsubspace}
\alias{plotsubspace}
\title{Plots covariance matrices}
\description{
  Represents covariance matrices as 3-d ellipsoids using the \code{rgl} package. 
Covariance matrices of dimension greater than 3 are plotted on the subspace 
defined by the first three eigenvectors.
}
\usage{
  plotsubspace(CA, CB=NULL, corr = FALSE, shadeCA = TRUE, 
      shadeCB = TRUE, axes.lab = FALSE, ...)
}
\arguments{
  \item{CA}{Matrix}
  \item{CB}{Optional second matrix}
  \item{corr}{If \code{TRUE} the covariance matrices are transformed into 
correlation matrices}
  \item{shadeCA}{If \code{TRUE} the ellipsoid is solid, if \code{FALSE} the 
ellipsoid is wireframe}
  \item{shadeCB}{If \code{TRUE} the ellipsoid is solid, if \code{FALSE} the 
ellipsoid is wireframe}
  \item{axes.lab}{If \code{TRUE} the axes are labelled with the eigenvectors}
  \item{...}{further arguments to be passed}
}
\details{
  The matrix CA is always red, and the matrix CB if given is always blue. The 
subspace is defined by the first three eigenvectors of CA, and the percentage of 
variance for each matrix along these three dimensions is given in the plot 
title. 
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk} with code taken from the \code{rgl} package}
\seealso{\code{\link[rgl]{rgl}}}
\examples{
 if(requireNamespace("rgl")!=FALSE){
   G1<-rIW(diag(4),10)
   G2<-G1*1.2
 #  plotsubspace(G1, G2, shadeCB=FALSE)
 # commented out because of problems with rgl 
 } 
}
\keyword{hplot}

