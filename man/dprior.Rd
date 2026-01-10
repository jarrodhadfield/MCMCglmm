\name{dprior}
\alias{dprior}
\title{Prior Density of Variance or Standard Deviation}
\description{Prior density of variance or standard deviation as specified in a MCMCglmm prior}
\usage{
dprior(x, prior, sd=FALSE, log=FALSE)
}
\arguments{
 \item{x}{vector of observations}
 \item{prior}{list: with elements \code{V}, \code{nu} and (optionally) \code{alpha.mu} and \code{alpha.V}}
 \item{sd}{logical: if TRUE the density of the standard deviation is returned}
 \item{log}{logical: if TRUE the log density is returned}
}
\value{
  numeric 
}
\details{

If \code{alpha.V} is zero (or \code{NULL}) the density for the variance is a scalar inverse-Wishart (i.e. inverse-gamma) distribution. If \code{alpha.V} is non-zero \code{MCMCglmm} uses parameter expansion and the density for the variance is a scaled non-central F distribution with the numerator degrees of freedom equal to one. The density for the standard deviation is a folded scaled non-central t distribution. 

}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{rIW}}
\examples{

prior<-list(V=1, nu=1, alpha.mu=0, alpha.V=10)
# parameter expanded prior

x<-seq(1e-6, 10, lenhth=100)

plot(dprior(x, prior)~x, type="l", xlab="Variance", ylab="Density")
# Density for variance: central F_{1,1} with scale 10

plot(dprior(sqrt(x), prior, sd=TRUE)~sqrt(x), type="l", xlab="Standard Deviation", ylab="Density")
# Density for standard deviation: half-Cauchy (t_1) with scale sqrt(10)

}
\keyword{distribution}
