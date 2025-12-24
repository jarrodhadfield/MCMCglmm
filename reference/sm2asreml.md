# Converts sparseMatrix to asreml's giv format

Converts sparseMatrix to asreml's giv format: row-ordered, upper
triangle sparse matrix.

## Usage

``` r
sm2asreml(A=NULL, rownames=NULL)
```

## Arguments

- A:

  sparseMatrix

- rownames:

  rownames of A

## Value

data.frame: if `A` was formed from a pedigree equivalent to giv format
returned by `asreml.Ainverse`

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

inverseA

## Examples

``` r
data(bird.families)
A<-inverseA(bird.families)
Aasreml<-sm2asreml(A$Ainv, A$node.names)
```
