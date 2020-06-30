\name{simulate.MCMCglmm}
\alias{simulate.MCMCglmm}
\title{Simulate method for GLMMs fitted with MCMCglmm}
\description{Simulated response vectors for GLMMs fitted with MCMCglmm}
\usage{
\method{simulate}{MCMCglmm}(object, nsim = 1, seed = NULL, newdata=NULL, marginal = object$Random$formula, 
          type = "response", it=NULL, posterior = "all", verbose=FALSE, \dots)
}
\arguments{
 \item{object}{an object of class \code{"MCMCglmm"}}
 \item{nsim}{number of response vectors to simulate.  Defaults to \code{1}.}
 \item{seed}{Either \code{NULL} or an integer that will be used in a call to \code{set.seed} before simulating the response vectors.  The default, \code{NULL} will not change the random generator state.}
 \item{newdata}{An optional data frame for which to simulate new observations}
 \item{marginal}{formula defining random effects to be maginalised}
 \item{type}{character; either "terms" (link scale) or "response" (data scale)}
 \item{it}{integer; optional, MCMC iteration on which predictions should be based}
 \item{posterior}{character; if \code{it} is  \code{NULL} should the response vector be simulated using the marginal posterior means ("mean") of the parameters, or the  posterior modes ("mode"), random draws from the posterior with replacement ("distribution") or  without replacement ("all")}
 \item{verbose}{logical;  if \code{TRUE}, warnings are issued with newdata when the original model has fixed effects that do not appear in newdata and/or newdata has random effects not present in the original model.}
 \item{\dots}{Further arguments to be passed}
}
\value{A matrix (with nsim columns) of simulated response vectors}  
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCglmm}}}
\keyword{models}


