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
#> [1] 0.8977743
#> 
#> $angles
#> [1] 19.19753 85.59436
#> 
#> $bisector
#>             [,1]        [,2]
#> [1,] -0.42593544 -0.05263835
#> [2,]  0.43626925  0.79678819
#> [3,] -0.78556362  0.51139583
#> [4,] -0.05219510 -0.05900718
#> [5,]  0.09172575  0.31201001
#> 
krzanowski.test(CA, CA, vecsA=1:2, vecsB=1:2)
#> $sumofS
#> [1] 2
#> 
#> $angles
#> [1] 0.000000e+00 1.478779e-06
#> 
#> $bisector
#>             [,1]        [,2]
#> [1,] -0.05102802 -0.59444276
#> [2,] -0.04330833  0.70394659
#> [3,] -0.95301180 -0.06227296
#> [4,] -0.23225098  0.35975604
#> [5,]  0.18261579  0.13339670
#> 
```
