\name{MCMCglmm-package}
\alias{MCMCglmm-package}
\docType{package}
\title{Multivariate Generalised Linear Mixed Models}
\description{
 MCMCglmm is a package for fitting Generalised Linear Mixed Models using Markov 
chain Monte Carlo techniques (Hadfield 2009). Most commonly used distributions like the normal 
and the Poisson are supported together with some useful but less popular ones 
like the zero-inflated Poisson and the multinomial. Missing values and left, 
right and interval censoring are accommodated for all traits. The package also 
supports multi-trait models where the multiple responses can follow different 
types of distribution.  The package allows various residual and random-effect 
variance structures to be specified including heterogeneous variances, 
unstructured covariance matrices and random regression (e.g. random slope 
models).  Three special types of variance structure that can be specified are 
those associated with pedigrees (animal models), phylogenies (the comparative 
method) and measurement error (meta-analysis).

The package makes heavy use of results in Sorensen & Gianola (2002) and Davis 
(2006) which taken together result in what is hopefully a fast and efficient 
routine. Most small to medium sized problems should take seconds to a few 
minutes, but large problems (> 20,000 records) are possible.  My interest is in 
evolutionary biology so there are also several functions for applying Rice's 
(2004) tensor analysis to real data and functions for visualising and comparing 
matrices.

Please read the tutorial \code{vignette("Tutorial", "MCMCglmm")} or the course notes
 \code{vignette("CourseNotes", "MCMCglmm")} 
}
\references{
	Hadfield, J.D. (2009) MCMC methods for Multi-response Generalised 
Linear Mixed Models: The MCMCglmm R Package, \emph{submitted}

	Sorensen, D. & Gianola, D. (2002) Likelihood, Bayesian and MCMC Methods in Quantitative 
Genetics, Springer

	Davis, T.A. (2006) Direct Methods for Sparse Linear Systems, SIAM
	
	Rice (2004) Evolutionary Theory: Mathematical and Conceptual Foundations, Sinauer
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}


