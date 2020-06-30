\name{summary.MCMCglmm}
\alias{summary.MCMCglmm}
\alias{print.summary.MCMCglmm}
\alias{print.MCMCglmm}
\title{Summarising GLMM Fits from MCMCglmm}
\description{\code{summary} method for class \code{"MCMCglmm"}. The returned object is suitable for printing with the \code{print.summary.MCMCglmm} method.
}
\usage{
\method{summary}{MCMCglmm}(object, random=FALSE, \dots)
}
\arguments{
 \item{object}{an object of class \code{"MCMCglmm"}}
 \item{random}{logical: should the random effects be summarised}
 \item{\dots}{Further arguments to be passed}
}
\value{
  \item{DIC}{Deviance Information Criterion}  
  \item{fixed.formula}{model formula for the fixed terms}
  \item{random.formula}{model formula for the random terms}
  \item{residual.formula}{model formula for the residual terms}
  \item{solutions}{posterior mean, 95\% HPD interval, MCMC p-values and effective sample size of fixed (and random) effects}
  \item{Gcovariances}{posterior mean, 95\% HPD interval and effective sample size of random effect (co)variance components}
  \item{Gterms}{indexes random effect (co)variances by the component terms defined in the random formula}
  \item{Rcovariances}{posterior mean, 95\% HPD interval and effective sample size of residual (co)variance components}
  \item{Rterms}{indexes residuals (co)variances by the component terms defined in the rcov formula}
  \item{csats}{chain length, burn-in and thinning interval}
  \item{cutpoints}{posterior mean, 95\% HPD interval and effective sample size of cut-points from an ordinal model}
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCglmm}}}
\keyword{models}

