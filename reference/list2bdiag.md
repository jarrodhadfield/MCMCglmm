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
#>           1         2         3         1         2
#> 1 1.2890318 0.6395845 0.4618807 0.0000000 0.0000000
#> 2 0.6395845 2.0298096 1.1601762 0.0000000 0.0000000
#> 3 0.4618807 1.1601762 2.5573003 0.0000000 0.0000000
#> 1 0.0000000 0.0000000 0.0000000 0.9667018 0.1474586
#> 2 0.0000000 0.0000000 0.0000000 0.1474586 0.9796461
```
