\name{rIW}
\alias{rIW}
\title{Random Generation from the Conditional Inverse Wishart Distribution}
\description{Samples from the inverse Wishart distribution, with the possibility 
of conditioning on a diagonal submatrix}
\usage{
rIW(V, nu, fix=NULL, n=1, CM=NULL)
}
\arguments{
 \item{V}{Expected (co)varaince matrix as \code{nu} tends to infinity}
 \item{nu}{degrees of freedom}
 \item{fix}{optional integer indexing the partition to be conditioned on}
 \item{n}{integer: number of samples to be drawn}
 \item{CM}{matrix: optional matrix to condition on. If not given, and \code{fix!=NULL},  V_22 is conditioned on}
}
\details{
If \eqn{{\bf W^{-1}}}{solve(W)} is a draw from the inverse Wishart,  \code{fix} indexes the diagonal element of \eqn{{\bf W^{-1}}}{solve(W)}  which partitions \eqn{{\bf W^{-1}}}{solve(W)} into 4 submatrices. \code{fix} indexes the upper left corner of the lower 
diagonal matrix and it is this matrix that is conditioned on.
 
 For example partioning \eqn{{\bf W^{-1}}}{solve(W)} such that 

 \deqn{
{\bf W^{-1}} = \left[
             \begin{array}{cc}
		      {\bf W^{-1}}_{11}&{\bf W^{-1}}_{12}\\
			   {\bf W^{-1}}_{21}&{\bf W^{-1}}_{22}\\
             \end{array}
		   \right]
  }{solve(W) = solve(W)_11 solve(W)_12}
  \deqn{}{          solve(W)_21 solve(W)_22}

 fix indexes the upper left corner of \eqn{{\bf W^{-1}}_{22}}{solve(W)_22}. If \code{CM!=NULL} then \eqn{{\bf W^{-1}}_{22}}{solve(W)_22} is fixed at \code{CM}, otherwise \eqn{{\bf W^{-1}}_{22}}{solve(W)_22} is fixed at  \eqn{\texttt{V}_{22}}{V_22}. For example, if \code{dim(V)}=4 and \code{fix=2} then \eqn{{\bf W^{-1}}_{11}}{solve(W)_11} is a 1X1 matrix and \eqn{{\bf W^{-1}}_{22}}{solve(W)_22}  is a 3X3 matrix.
}
\value{
  if \code{n} = 1 a matrix equal in dimension to \code{V}, if \code{n}>1 a 
matrix of dimension \code{n} x \code{length(V)}
}
\note{In versions of MCMCglmm >1.10 the arguments to \code{rIW} have changed so that they are more intuitive in the context of \code{\link{MCMCglmm}}.  Following the notation of Wikipedia (\url{https://en.wikipedia.org/wiki/Inverse-Wishart_distribution}) the inverse scale matrix \eqn{{\bm \Psi}=(\texttt{V*nu})}{Psi = \code{V*nu}}. In earlier versions of MCMCglmm (<1.11) \eqn{{\bm \Psi} = \texttt{V}^{-1}}{Psi=\code{solve(V)}}. Although the old parameterisation is consistent with the \code{\link[MCMCpack]{riwish}} function in MCMCpack and the \code{\link[bayesm]{rwishart}} function in bayesm it is inconsistent with the prior definition for \code{\link{MCMCglmm}}.  The following pieces of code are sampling from the same distributions: 

\tabular{ll}{
\code{\link[MCMCpack]{riwish}(nu, nu*V)} \tab from MCMCpack\cr
\code{\link[bayesm]{rwishart}(nu, solve(nu*V))$IW} \tab from bayesm\cr
\code{rIW(nu, solve(nu*V))} \tab from MCMCglmm <1.11\cr
\code{rIW(V, nu)} \tab from MCMCglmm >=1.11\cr
}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\references{Korsgaard, I.R. et. al. 1999 Genetics Selection Evolution 31 (2) 
177:181}
\seealso{\code{\link[bayesm]{rwishart}}, \code{\link[MCMCpack]{rwish}}}
\examples{
nu<-10
V<-diag(4)
rIW(V, nu, fix=2)
}
\keyword{distribution}

