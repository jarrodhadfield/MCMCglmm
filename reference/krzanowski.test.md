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
#> [1] 1.210027
#> 
#> $angles
#> [1] 19.16199 55.68724
#> 
#> $bisector
#>             [,1]          [,2]
#> [1,] -0.29956438 -0.2788995867
#> [2,]  0.55409313  0.3058949472
#> [3,] -0.67188154  0.7026849270
#> [4,] -0.37870132 -0.5786852372
#> [5,]  0.09166514 -0.0007690872
#> 
krzanowski.test(CA, CA, vecsA=1:2, vecsB=1:2)
#> Warning: NaNs produced
#> $sumofS
#> [1] 2
#> 
#> $angles
#> [1]          NaN 1.478779e-06
#> 
#> $bisector
#>             [,1]        [,2]
#> [1,] -0.49671981 -0.15079396
#> [2,]  0.42146712 -0.44890151
#> [3,]  0.04445777  0.87512346
#> [4,] -0.65331641 -0.08603530
#> [5,]  0.38319196  0.05005473
#> 
```
