\name{resolve_prior}
\alias{resolve_prior}
\title{Resolves a Prior Specification in MCMCglmm}
\description{Generates the (co)variance matrix prior parameters in MCMCglmm given a prior specification that might utilise a \code{prior_generator}.}
\usage{
  resolve_prior(prior_element, k=NULL)
}
\arguments{
 \item{prior_element}{prior spefictaion for a (co)variance matrix}
 \item{k}{dimension of the (co)variance matrix}
 \item{vtype}{type of covariance structure used: currently only \code{us} and \code{idh} are allowed}
}
\details{
  Primarily used internally to generate the (co)variance matrix prior parameters used in MCMCglmm: \code{V}, \code{nu}, \code{alpha.mu} and \code{alpha.V}. If these are already passed as a list, the list is simply returned. If a \code{prior_generator} is used then the prior parameters are generated given dimension of the (co)variance matrix, \code{k} and its stucture type \code{vtype}.
}
\value{
  list of prior parameters for MCMCglmm
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{prior_generator}},\code{\link{IW}}, \code{\link{IG}}, \code{\link{F}}, \code{\link{tSD}}}

\examples{
prior_element1<-list(V=diag(2)/2, nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*1000)

resolve_prior(prior_element1)

prior_element2<-F(df2=1, scale=1000)

resolve_prior(prior_element2, k=2, vtype="us")

resolve_prior(prior_element2, k=2, vtype="idh")

}
\keyword{distribution}

