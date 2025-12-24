# Krzanowski's Comparison of Subspaces

Calculates statistics of Krzanowski's comparison of subspaces.

## Usage

``` r
krzanowski.test(CA, CB, vecsA, vecsB, corr = FALSE, ...)
```

## Arguments

- CA:

  Matrix A

- CB:

  Matrix B

- vecsA:

  Vector of integers indexing the eigenvectors determining the subspace
  of A

- vecsB:

  Vector of integers indexing the eigenvectors determining the subspace
  of B

- corr:

  logical; if `TRUE` the variances of A and B are standardised

- ...:

  further arguments to be passed

## Value

- sumofS:

  metric for overall similarity with 0 indicting no similarity and a
  value of `length(vecsA)` for identical subspaces

- angles:

  angle in degrees between each best matched pair of vectors

- bisector:

  vector that lies between each best matched pair of vectors

## References

Krzanowski, W.J. (2000) Principles of Multivariate Analysis. OUP

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
CA<-rIW(diag(5),10, n=1)
CB<-rIW(diag(5),10, n=1)
krzanowski.test(CA, CB, vecsA=1:2, vecsB=1:2)
#> $sumofS
#> [1] 0.4755277
#> 
#> $angles
#> [1] 48.66022 78.57510
#> 
#> $bisector
#>             [,1]       [,2]
#> [1,]  0.61286510  0.6305307
#> [2,]  0.61437157 -0.5589670
#> [3,] -0.35097369 -0.2602179
#> [4,]  0.35011737 -0.3530748
#> [5,] -0.03434013  0.3124287
#> 
krzanowski.test(CA, CA, vecsA=1:2, vecsB=1:2)
#> $sumofS
#> [1] 2
#> 
#> $angles
#> [1] 0.000000e+00 1.478779e-06
#> 
#> $bisector
#>            [,1]        [,2]
#> [1,] -0.5224690 -0.05079276
#> [2,] -0.2890904  0.41888305
#> [3,]  0.7625666 -0.04001207
#> [4,]  0.1082726  0.89711567
#> [5,] -0.2241029 -0.12465791
#> 
```
