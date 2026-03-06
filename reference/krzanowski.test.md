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
#> [1] 0.9271759
#> 
#> $angles
#> [1] 27.20786 68.34061
#> 
#> $bisector
#>            [,1]        [,2]
#> [1,] -0.1983429  0.39576182
#> [2,] -0.1027287  0.87848379
#> [3,]  0.4895669 -0.08465776
#> [4,] -0.8353376 -0.25352035
#> [5,] -0.1124375  0.01411782
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
#> [1,] -0.1263269 -0.16113845
#> [2,] -0.2089699  0.96731263
#> [3,]  0.6044028  0.08468567
#> [4,] -0.7553711 -0.16811378
#> [5,]  0.0669688  0.05391444
#> 
```
