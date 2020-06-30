\name{me}
\alias{me}
\title{Design Matrix for Measurement Error Model}
\description{Sets up design matrix for measurement error models.}
\usage{
me(formula, error=NULL, group=NULL, type="classical")

}
\arguments{
 \item{formula}{\code{\link{formula}} for the fixed effects.}
 \item{error}{character; name of column in \code{\link{data.frame}} in which standard error (\code{type="classical"} or \code{type="berkson"}) or miscalssification error  (\code{type="dclassical"}) is stored.}
 \item{group}{name of column in \code{\link{data.frame}} in which groups are stored. Rows of the design matrix with the same \code{group} level are assumed to pertain to the same obsevation of the covariate that is measured with error.}
  \item{type}{character; one of \code{type="classical"}, \code{type="berkson"}, \code{type="dclassical"} or \code{type="dberkson"} (see details)}
 }
\value{
  design matrix, with a prior distribution attribute 
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\keyword{models}

