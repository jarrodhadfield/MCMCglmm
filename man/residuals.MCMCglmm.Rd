\name{residuals.MCMCglmm}
\alias{residuals.MCMCglmm}
\title{Residuals form a GLMM fitted with MCMCglmm}
\description{\code{residuals} method for class \code{"MCMCglmm"}.
}
\usage{
\method{residuals}{MCMCglmm}(object, type = c("deviance", "pearson", "working",
                                "response", "partial"), \dots)
}
\arguments{
 \item{object}{an object of class \code{"MCMCglmm"}}
 \item{type}{ the type of residuals which should be returned. The alternatives are: \code{"deviance"} (default), \code{"pearson"},\code{"working"}, \code{"response"}, and \code{"partial"}.}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of residuals
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{residuals}}, \code{\link{MCMCglmm}}}
\keyword{models}

