\name{spl}
\alias{spl}
\title{Orthogonal Spline Design Matrix}
\description{Orthogonal Spline Design Matrix}
\usage{
spl(x,  k=10, knots=NULL, type="LRTP")
}
\arguments{
  \item{x}{a numeric covariate}
  \item{k}{integer, defines knot points at the 1:k/(k+1) quantiles of x}
  \item{knots}{vector of knot points}
  \item{type}{type of spline - currently only low-rank thin-plate ("LRTP") are implemented}
}
\value{
	Design matrix post-multiplied by the inverse square root of the penalty matrix 
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
\dontrun{
x<-rnorm(100)
y<-x^2+cos(x)-x+0.2*x^3+rnorm(100)
plot(y~x)
lines((x^2+cos(x)-x+0.2*x^3)[order(x)]~sort(x))

dat<-data.frame(y=y, x=x)

m1<-MCMCglmm(y~x, random=~idv(spl(x)), data=dat, pr=TRUE, verbose=FALSE) # penalised smoother
m2<-MCMCglmm(y~x+spl(x),data=dat,  verbose=FALSE)                        # non-penalised

pred1<-(cbind(m1$X,m1$Z)\%*\%colMeans(m1$Sol))@x
pred2<-(cbind(m2$X)\%*\%colMeans(m2$Sol))@x

lines(pred1[order(x)]~sort(x), col="red")
lines(pred2[order(x)]~sort(x), col="green")

m1$DIC-mean(m1$Deviance)  # effective number of parameters < 13
m2$DIC-mean(m2$Deviance)  # effective number of parameters ~ 13
}
}
\keyword{manip}

