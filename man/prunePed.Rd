\name{prunePed}
\alias{prunePed}
\title{Pedigree pruning}
\description{
 Creates a subset of a pedigree by retaining the ancestors of a specified subset 
of individuals
}
\usage{
prunePed(pedigree, keep, make.base=FALSE)
}
\arguments{
  \item{pedigree}{pedigree with id in column 1 dam in column 2 and sire in 
column 3}
  \item{keep}{individuals in pedigree for which the ancestors should be 
retained}
  \item{make.base}{logical: should ancestors that do not provide additional information be discarded?}  
}
\value{
  subsetted pedigree
}
\note{
   If the individuals in \code{keep} are the only phenotyped individuals for some analysis then some non-phenotyped individuals can often be discarded if they are not responsible for pedigree links between phenotyped individuals.  In the simplest case (\code{make.base=FALSE}) all ancestors of phenotyped individuals will be retained, although further pruning may be possible using \code{make.base=TRUE}. In this case all pedigree links that do not connect phenotyped individuals are discarded resulting in some individuals becoming part of the base population. In terms of variance component and fixed effect estimation pruning the pedigree should have no impact on the target posterior distribution, although convergence and mixing may be better because there is less missing data. 
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk} + Michael Morrissey}
\keyword{manip}


