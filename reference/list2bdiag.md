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
#>             1           2          3          1          2
#> 1  1.69461353 -0.07141992  0.9631322  0.0000000  0.0000000
#> 2 -0.07141992  0.99251923 -0.4896201  0.0000000  0.0000000
#> 3  0.96313221 -0.48962009  1.5795523  0.0000000  0.0000000
#> 1  0.00000000  0.00000000  0.0000000  4.2599594 -0.3691509
#> 2  0.00000000  0.00000000  0.0000000 -0.3691509  0.6253375
```
