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
#> [1] 1.070153
#> 
#> $angles
#> [1] 27.21882 58.09290
#> 
#> $bisector
#>             [,1]       [,2]
#> [1,] -0.09649454 -0.1965638
#> [2,]  0.14095438 -0.1391929
#> [3,] -0.08629227 -0.5994571
#> [4,] -0.24268854 -0.7250194
#> [5,]  0.95103974 -0.2387176
#> 
krzanowski.test(CA, CA, vecsA=1:2, vecsB=1:2)
#> Warning: NaNs produced
#> $sumofS
#> [1] 2
#> 
#> $angles
#> [1]          NaN 2.957559e-06
#> 
#> $bisector
#>            [,1]        [,2]
#> [1,] -0.5536261  0.05845951
#> [2,]  0.1037437 -0.29794082
#> [3,] -0.3389537  0.20448050
#> [4,] -0.7189961  0.10568850
#> [5,] -0.2255892 -0.92457094
#> 
```
