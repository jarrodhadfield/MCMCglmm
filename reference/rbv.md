# Random Generation of MVN Breeding Values and Phylogenetic Effects

Random Generation of MVN Breeding Values and Phylogenetic Effects

## Usage

``` r
rbv(pedigree, G, nodes="ALL", scale=TRUE, ggroups=NULL, gmeans=NULL)
```

## Arguments

- pedigree:

  ordered pedigree with 3 columns id, dam and sire or a `phylo` object.

- G:

  (co)variance matrix

- nodes:

  effects for pedigree/phylogeny nodes to be returned. The default,
  `nodes="ALL"` returns effects for all individuals in a pedigree or
  nodes in a phylogeny (including ancestral nodes). For phylogenies
  `nodes="TIPS"` returns effects for the tips only, and for pedigrees a
  vector of ids can be passed to `nodes` specifying the subset of
  individuals for which animal effects are returned.

- scale:

  logical: should a phylogeny (needs to be ultrametric) be scaled to
  unit length (distance from root to tip)?

- ggroups:

  optional; vector of genetic groups

- gmeans:

  matrix of mean breeding value for genetic groups (rows) by traits
  (columns)

## Value

matrix of breeding values/phylogenetic effects

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
data(bird.families)
bv<-rbv(bird.families, diag(2))
```
