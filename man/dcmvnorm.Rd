\name{dcmvnorm}
\alias{dcmvnorm}
\title{Density of a (conditional) multivariate normal variate}
\description{Density of a (conditional) multivariate normal variate}
\usage{
dcmvnorm(x, mean = 0, V = 1, keep=1, cond=(1:length(x))[-keep], log=FALSE)
}
\arguments{
 \item{x}{vector of observations}
 \item{mean}{vector of means}
 \item{V}{covariance matrix}
 \item{keep}{vector of integers: observations for which density is required}
 \item{cond}{vector of integers: observations to condition on}
 \item{log}{if TRUE, density p is given as log(p)}
}
\value{
  numeric 
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
V1<-cbind(c(1,0.5), c(0.5,1))
dcmvnorm(c(0,2), c(0,0), V=V1, keep=1, cond=2)
# density of x[1]=0 conditional on x[2]=2 given 
# x ~ MVN(c(0,0), V1) 

dcmvnorm(c(0,2), c(0,0), V=V1, keep=1, cond=NULL)
# density of x[1]=0 marginal to x[2] 
dnorm(0,0,1)
# same as univariate density 

V2<-diag(2)
dcmvnorm(c(0,2), c(0,0), V=V2, keep=1, cond=2)
# density of x[1]=0 conditional on x[2]=2 given 
# x ~ MVN(c(0,0), V2) 
dnorm(0,0,1)
# same as univariate density because V2 is diagonal
}
\keyword{distribution}
