# Pedigree pruning

Creates a subset of a pedigree by retaining the ancestors of a specified
subset of individuals

## Usage

``` r
prunePed(pedigree, keep, make.base=FALSE)
```

## Arguments

- pedigree:

  pedigree with id in column 1 dam in column 2 and sire in column 3

- keep:

  individuals in pedigree for which the ancestors should be retained

- make.base:

  logical: should ancestors that do not provide additional information
  be discarded?

## Value

subsetted pedigree

## Note

If the individuals in `keep` are the only phenotyped individuals for
some analysis then some non-phenotyped individuals can often be
discarded if they are not responsible for pedigree links between
phenotyped individuals. In the simplest case (`make.base=FALSE`) all
ancestors of phenotyped individuals will be retained, although further
pruning may be possible using `make.base=TRUE`. In this case all
pedigree links that do not connect phenotyped individuals are discarded
resulting in some individuals becoming part of the base population. In
terms of variance component and fixed effect estimation pruning the
pedigree should have no impact on the target posterior distribution,
although convergence and mixing may be better because there is less
missing data.

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk> + Michael Morrissey
