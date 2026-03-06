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
#>             1           2          3            1            2
#> 1  0.81643550  0.08136069 -0.4021869  0.000000000  0.000000000
#> 2  0.08136069  0.95855318 -0.2277594  0.000000000  0.000000000
#> 3 -0.40218690 -0.22775945  2.0679720  0.000000000  0.000000000
#> 1  0.00000000  0.00000000  0.0000000  0.939165765 -0.008854768
#> 2  0.00000000  0.00000000  0.0000000 -0.008854768  0.823375164
```
