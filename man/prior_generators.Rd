\name{prior_generators}
\alias{prior_generators}
\alias{IW}
\alias{IG}
\alias{F}
\alias{tSD}
\title{Generator Functions for Priors in MCMCglmm}
\description{Functions for generating (co)variance matrix prior specifications in MCMCglmm that result in specified inverse-Wishart, inverse-gamma or central-\eqn{F}{F} marginal priors for the variances.}
\usage{
  IW(V=1, nu=0.002)

  IG(shape=0.001, scale=0.001)
  
  F(df2=1, scale=1000)
  
  tSD(df=1, scale=sqrt(1000))
}
\arguments{
 \item{V}{expected varaince as \code{nu} tends to infinity in a scalar inverse-Wishart prior}
 \item{nu}{degrees of freedom in a scalar inverse-Wishart prior}
 \item{shape}{shape parameter of the inverse-gamma prior}
 \item{scale}{scale parameter of the inverse-gamma, the scaled-F or scaled half-t}
 \item{df}{degrees of freedom for the half-t prior on the standard deviation}
 \item{df2}{denominator degrees of freedom for F prior (numerator degree-of-freedom is one)}
}
\details{
  Each genertor function returns a function that generates a list of prior arguments need to specific (co)variance matrix priors in \code{MCMCglmm}. Those prior arguments result in the marginal distributions for the variancess  being those specified in generator function. Since the appropriate prior arguments depend on the dimension of the (co)variance matrix, they are evalauted at run time once the dimension is determined using \code{resove_prior}.  
}
\value{
  function of class \code{prior_generator}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{prior_generator}}}
\examples{
resolve_prior(F(df2=1, scale=1000), k=2)
}
\keyword{distribution}

