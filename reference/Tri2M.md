# Lower/Upper Triangle Elements of a Matrix

Lower/Upper triangle elements of a matrix or forms a matrix from a
vector of lower/upper triangle elements

## Usage

``` r
Tri2M(x, lower.tri = TRUE, reverse = TRUE, diag = TRUE)
```

## Arguments

- x:

  Matrix or vector

- lower.tri:

  If `x` is a matrix then the lower triangle (`TRUE`) or upper triangle
  `FALSE` elements (including diagonal elements) are returned. If `x` is
  a vector a matrix is formed under the assumption that `x` are the
  lower triangle (`TRUE`) or upper triangle (`FALSE`) elements.

- reverse:

  logical: if`TRUE` a symmetric matrix is formed, if `FALSE` the
  remaining triangle is left as zeros.

- diag:

  logical: if`TRUE` diagonal elements are included.

## Value

numeric or matrix

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
M<-rIW(diag(3), 10)
x<-Tri2M(M)
x
#> [1]  1.9721942 -0.2289415  0.2278035  1.1638801 -0.1994067  0.8686290
Tri2M(x, reverse=TRUE)
#>            [,1]       [,2]       [,3]
#> [1,]  1.9721942 -0.2289415  0.2278035
#> [2,] -0.2289415  1.1638801 -0.1994067
#> [3,]  0.2278035 -0.1994067  0.8686290
Tri2M(x, reverse=FALSE)
#>            [,1]       [,2]     [,3]
#> [1,]  1.9721942  0.0000000 0.000000
#> [2,] -0.2289415  1.1638801 0.000000
#> [3,]  0.2278035 -0.1994067 0.868629
```
