\name{buildV}
\alias{buildV}
\title{Forms expected (co)variances for GLMMs fitted with MCMCglmm}
\description{Forms the expected covariance structure of link-scale observations for GLMMs fitted with MCMCglmm}
\usage{
buildV(object, marginal=object$Random$formula, diag=TRUE, it=NULL, posterior="mean", \dots)
}
\arguments{
 \item{object}{an object of class \code{"MCMCglmm"}}
 \item{marginal}{formula defining random effects to be maginalised}
 \item{diag}{logical; if \code{TRUE} the covariances betwween observations are not calculated}
 \item{it}{integer; optional, MCMC iteration on which covariance structure should be based}
 \item{posterior}{character; if \code{it} is  \code{NULL} should the covariance structure be based on the marginal posterior means (\code{'mean'}) of the VCV parameters, or the  posterior modes (\code{'mode'}), or a random draw from the posterior with replacement (\code{'distribution'}).  If \code{posterior=="all"} the posterior distribution of observation variances is returned}
 \item{\dots}{Further arguments to be passed}
}
\value{If \code{diag=TRUE} an n by n covariance matrix. If \code{diag=FALSE} and  \code{posterior!="all"} an 1 by n matrix of variances. If \code{posterior=="all"} an nit by n matrix of variances (where nit is the number of saved MCMC iterations).}  
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCglmm}}}
\keyword{models}


