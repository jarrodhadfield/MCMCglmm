\name{rprior}
\alias{rprior}
\title{Simulate Covariance Matrices from Prior}
\description{Simulates covariance matrices from prior as specified in a MCMCglmm prior}
\usage{
rprior(n, prior, vtype="us")
}
\arguments{
 \item{n}{number of observations}
 \item{prior}{list: with elements \code{V}, \code{nu} and (optionally) \code{alpha.mu} and \code{alpha.V}}
 \item{vtype}{character: variance structure type with default \code{us}}
}
\value{
  numeric 
}
\details{

If \code{alpha.V} is \code{NULL} the distribution is an inverse-Wishart distribution (i.e. inverse-gamma in the univariate case). If \code{alpha.V} is non-null, \code{MCMCglmm} uses parameter expansion and the distribution is an inverse-Wishart mixture with no closed form density function expect in the univariate case (scaled F with 1 numerator degree of freedom). 

}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{rIW}, \code{rprior}}
\examples{

prior<-resolve_prior(F(10, 20), k=3, vtype="us")
# parameter expanded prior for 3x3 covariance matrix with scaled (20) central F_{1,10} marginal variances

V<-rprior(2000, prior)

hist(V[,1], freq=FALSE, breaks=100, main="", xlab="Variance")
x<-seq(1e-6, max(V[,1]), length=1000)

mprior<-resolve_prior(F(10, 20), k=1, vtype="us")
# univariate prior for central F_{1,10} with scale 20

lines(dprior(x, mprior)~x)
# Density for variance


}
\keyword{distribution}
