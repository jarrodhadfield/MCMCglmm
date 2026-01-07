# Forms the direct sum from a list of matrices

Forms a block-diagonal matrix from a list of matrices

## Usage

``` r
list2bdiag(x)
```

## Arguments

- x:

  list of square matrices

## Value

matrix

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
M<-list(rIW(diag(3), 10), rIW(diag(2), 10))
list2bdiag(M)
#>            1          2          3            1            2
#> 1  0.9125119  0.1479412 -0.1243541  0.000000000  0.000000000
#> 2  0.1479412  0.5327244 -0.1697929  0.000000000  0.000000000
#> 3 -0.1243541 -0.1697929  2.0679720  0.000000000  0.000000000
#> 1  0.0000000  0.0000000  0.0000000  0.939165765 -0.008854768
#> 2  0.0000000  0.0000000  0.0000000 -0.008854768  0.823375164
```
