\name{path}
\alias{path}
\title{Design Matrix for Path Analyses}
\description{Forms design matrix for path analyses that involve paths within residual blocks}
\usage{
path(cause=NULL, effect=NULL, k)
}
\arguments{
  \item{cause}{integer; index of predictor `trait' within residual block}
  \item{effect}{integer; index of response `trait' within residual block }
  \item{k}{integer; dimension of residual block }
}
\value{
    design matrix 
}
\note{For more general path anlaytic models see \link{sir} which allows paths to exist between responses that are not in the same residual block. However, \code{sir} does not handle non-Gaussian or missing responses. Note that path models involving non-Gaussian data are defined on the link scale which may not always be appropriate.
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\link{sir}}
\examples{
path(1, 2,2)
}
\keyword{manip}

