\name{rbv}
\alias{rbv}
\title{Random Generation of MVN Breeding Values and Phylogenetic Effects}
\description{Random Generation of MVN Breeding Values and Phylogenetic Effects}
\usage{
rbv(pedigree, G, nodes="ALL", scale=TRUE, ggroups=NULL, gmeans=NULL)
}
\arguments{
 \item{pedigree}{ordered pedigree with 3 columns id, dam and sire or a 
\code{phylo} object.}
 \item{G}{(co)variance matrix}
  \item{nodes}{effects for pedigree/phylogeny nodes to be returned. The default, 
\code{nodes="ALL"} returns effects for all individuals in a pedigree or nodes in 
a phylogeny (including ancestral nodes). For phylogenies \code{nodes="TIPS"} 
returns effects for the tips only, and for pedigrees a vector of ids can be 
passed to \code{nodes} specifying the subset of individuals for which animal 
effects are returned.}
 \item{scale}{logical: should a phylogeny (needs to be ultrametric) be scaled to 
unit length (distance from root to tip)?}
 \item{ggroups}{optional; vector of genetic groups}
 \item{gmeans}{matrix of mean breeding value for genetic groups (rows) by traits (columns)}
}
\value{
  matrix of breeding values/phylogenetic effects
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
data(bird.families)
bv<-rbv(bird.families, diag(2))
}
\keyword{distribution}

