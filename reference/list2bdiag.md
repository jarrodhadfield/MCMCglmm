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
#>            1          2          3           1           2
#> 1 0.71371267 0.23201910 0.02017168 0.000000000 0.000000000
#> 2 0.23201910 0.97594747 0.07636249 0.000000000 0.000000000
#> 3 0.02017168 0.07636249 0.94680596 0.000000000 0.000000000
#> 1 0.00000000 0.00000000 0.00000000 3.067658855 0.005098866
#> 2 0.00000000 0.00000000 0.00000000 0.005098866 0.742262001
```
