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
#>           1         2         3           1           2
#> 1 1.1737345 0.4860404 0.1582004 0.000000000 0.000000000
#> 2 0.4860404 1.0110503 0.2553205 0.000000000 0.000000000
#> 3 0.1582004 0.2553205 1.0254333 0.000000000 0.000000000
#> 1 0.0000000 0.0000000 0.0000000 1.252368336 0.004165354
#> 2 0.0000000 0.0000000 0.0000000 0.004165354 0.418904231
```
