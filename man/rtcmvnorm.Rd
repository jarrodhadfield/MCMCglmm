\name{rtcmvnorm}
\alias{rtcmvnorm}
\title{Random Generation from a Truncated Conditional Normal Distribution}
\description{Samples from the Truncated Conditional Normal Distribution}
\usage{
rtcmvnorm(n = 1, mean = 0, V = 1, x=0, keep=1, lower = -Inf, upper = Inf)
}
\arguments{
 \item{n}{integer: number of samples to be drawn}
 \item{mean}{vector of means}
 \item{V}{covariance matrix}
 \item{x}{vector of observations to condition on}
 \item{keep}{element of x to be sampled}
 \item{lower}{left truncation point}
 \item{upper}{right truncation point}
}
\value{
  vector 
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
par(mfrow=c(2,1))
V1<-cbind(c(1,0.5), c(0.5,1))
x1<-rtcmvnorm(10000, c(0,0), V=V1, c(0,2), keep=1, lower=-1, upper=1)
x2<-rtnorm(10000, 0, 1, lower=-1, upper=1)
plot(density(x1), main="Correlated conditioning observation")
lines(density(x2), col="red")
# denisties of conditional (black) and unconditional (red) distribution
# when the two variables are correlated (r=0.5) 

V2<-diag(2)
x3<-rtcmvnorm(10000, c(0,0), V=V2, c(0,2), keep=1, lower=-1, upper=1)
x4<-rtnorm(10000, 0, 1, lower=-1, upper=1)
plot(density(x3), main="Uncorrelated conditioning observation")
lines(density(x4), col="red")
# denisties of conditional (black) and unconditional (red) distribution
# when the two variables are uncorrelated (r=0) 

}
\keyword{distribution}

