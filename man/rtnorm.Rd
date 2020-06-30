\name{rtnorm}
\alias{rtnorm}
\title{Random Generation from a Truncated Normal Distribution}
\description{Samples from the Truncated Normal Distribution}
\usage{
rtnorm(n = 1, mean = 0, sd = 1, lower = -Inf, upper = Inf)
}
\arguments{
 \item{n}{integer: number of samples to be drawn}
 \item{mean}{vector of means}
 \item{sd}{vector of standard deviations}
 \item{lower}{left truncation point}
 \item{upper}{right truncation point}
}
\value{
  vector 
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\references{Robert, C.P. (1995) Statistics & Computing 5 121-125}
\seealso{\code{\link[msm]{rtnorm}}}
\examples{
hist(rtnorm(100, lower=-1, upper=1))
}
\keyword{distribution}

