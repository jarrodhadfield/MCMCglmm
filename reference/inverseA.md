# Inverse Relatedness Matrix and Phylogenetic Covariance Matrix

Henderson (1976) and Meuwissen and Luo (1992) algorithm for inverting
relatedness matrices, and Hadfield and Nakagawa (2010) algorithm for
inverting phylogenetic covariance matrices.

## Usage

``` r
inverseA(pedigree=NULL, nodes="ALL", scale=TRUE, reduced=FALSE,
     tol = .Machine$double.eps^0.5)
```

## Arguments

- pedigree:

  ordered pedigree with 3 columns: id, dam and sire, or a `phylo`
  object.

- nodes:

  `"ALL"` calculates the inverse for all individuals/nodes. For
  phylogenies `"TIPS"` calculates the inverse for the species tips only,
  and for pedigrees a vector of id's can be passed which inverts the
  relatedness matrix for that subset.

- scale:

  logical: should a phylogeny (needs to be ultrametric) be scaled to
  unit length (distance from root to tip)?

- reduced:

  logical: should childless nodes be dropped from the inverse and the
  pedigree/phylogeny representation be reduced?

- tol:

  numeric: differences in branch length smaller than this are ignored
  when assessing whether a tree is ultrametric.

## Value

- Ainv:

  inverse as `sparseMatrix`

- inbreeding:

  inbreeding coefficients/branch lengths

- pedigree:

  pedigree/pedigree representation of phylogeny

## References

Henderson, C.R. (1976) Biometrics 32 (1) 69:83

Quaas, R. L. and Pollak, E. J. (1980) Journal of Animal Science
51:1277-1287.

Meuwissen, T.H.E and Luo, Z. (1992) Genetic Selection Evolution 24 (4)
305:313

Hadfield, J.D. and Nakagawa, S. (2010) Journal of Evolutionary Biology
23 494-508

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
data(bird.families)
Ainv<-inverseA(bird.families)
```
