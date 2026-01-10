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
#> [1] 0.8280372
#> 
#> $angles
#> [1] 24.73400 86.80784
#> 
#> $bisector
#>           [,1]        [,2]
#> [1,] 0.3784868 -0.05115748
#> [2,] 0.6139955 -0.19118687
#> [3,] 0.1178499  0.86519304
#> [4,] 0.6742648  0.12146010
#> [5,] 0.1059981 -0.44443104
#> 
krzanowski.test(CA, CA, vecsA=1:2, vecsB=1:2)
#> $sumofS
#> [1] 2
#> 
#> $angles
#> [1] 0.000000e+00 1.909096e-06
#> 
#> $bisector
#>            [,1]        [,2]
#> [1,]  0.1340174 -0.50717840
#> [2,] -0.6610381 -0.59918400
#> [3,]  0.5221505 -0.05239904
#> [4,]  0.4733443 -0.61581706
#> [5,] -0.2199363  0.04209868
#> 
```
