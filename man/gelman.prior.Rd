\name{gelman.prior}
\alias{gelman.prior}
\title{Prior Covariance Matrix for Fixed Effects.}
\description{Prior covariance matrix for fixed effects had inputs been standardised as suggested in Gelman et al. (2008).}
\usage{
gelman.prior(formula, data, coef.scale=1, intercept.scale=coef.scale, singular.ok=FALSE)
}
\arguments{
 \item{formula}{\code{\link{formula}} for the fixed effects.}
 \item{data}{\code{\link{data.frame}}.}
 \item{coef.scale}{prior standard deviation for regression parameters  (had inputs been standardised): default=1.}
 \item{intercept.scale}{prior standard deviation for the intercept (had inputs been standardised): default=\code{coef.scale}.}
 \item{singular.ok}{logical: if \code{FALSE} linear dependencies in the fixed effects are removed. if \code{TRUE} they are left in an estimated, although all information comes form the prior}
  \item{\dots}{Further arguments to be passed}
 }
\value{
  list with elements \code{mu} (prior mean) and \code{V} prior covariance covariance matrix
 }
\references{
        Gelman, A. et al. (2008) The Annals of Appled Statistics 2 4  1360-1383
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\details{
 Gelman et al. (2008) suggest that the input variables in logistic regression are standardised and that the associated regression parameters are assumed independent in the prior.  Gelman et al. (2008) recommend a scaled t-distribution with a single degree of freedom (scaled Cauchy) and a scale of 10 for the intercept and 2.5 for the regression parameters. If the degree of freedom is infinity (i.e. a normal distribution) then a prior covariance matrix \code{B$V} can be defined for the regression parameters without input standardisation that corresponds to a diagonal prior \eqn{{\bf D}}{D} for the regression parameters had the inputs been standardised.  The diagonal elements of \eqn{{\bf D}}{D} are set to \code{coef.scale^2} except the first which is set to \code{intercept.scale^2}. Depending on the link-function, the presence of random effects, and the strength of prior required, suitable values for \code{coef.scale} and \code{intercept.scale} may differ than those recommened by Gelman et al. (2008). For details see \url{https://jarrodhadfield.github.io/MCMCglmm/course-notes/glm.html#gelman-prior-sec}.

}
\keyword{distribution}

