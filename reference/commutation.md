# Commutation Matrix

Forms an mn x mn commutation matrix which transforms \\vec({\bf A})\\
into \\vec({\bf A}^{'})\\, where \\{\bf A}\\ is an m x n matrix

## Usage

``` r
commutation(m, n)
```

## Arguments

- m:

  integer; number of rows of A

- n:

  integer; number of columns of A

## Value

Commutation Matrix

## References

Magnus, J. R. & Neudecker, H. (1979) Annals of Statistics 7 (2) 381-394

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
commutation(2,2)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    0    0    1    0
#> [3,]    0    1    0    0
#> [4,]    0    0    0    1
```
