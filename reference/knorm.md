# (Mixed) Central Moments of a Multivariate Normal Distribution

Forms a tensor of (mixed) central moments of a multivariate normal
distribution

## Usage

``` r
knorm(V, k)
```

## Arguments

- V:

  (co)variance matrix

- k:

  kth central moment, must be even

## Value

tensor

## References

Schott, J.R.(2003) Journal of Multivariate Analysis 87 (1) 177-190

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`dnorm`](https://rdrr.io/r/stats/Normal.html)

## Examples

``` r
V<-diag(2)
knorm(V,2)
#>       I2
#> I1     [,1] [,2]
#>   [1,]    1    0
#>   [2,]    0    1
#> attr(,"class")
#> [1] "tensor" "matrix"
knorm(V,4)
#> , , 1, 1
#> 
#>       I2
#> I1     [,1] [,2]
#>   [1,]    3    0
#>   [2,]    0    1
#> 
#> , , 2, 1
#> 
#>       I2
#> I1     [,1] [,2]
#>   [1,]    0    1
#>   [2,]    1    0
#> 
#> , , 1, 2
#> 
#>       I2
#> I1     [,1] [,2]
#>   [1,]    0    1
#>   [2,]    1    0
#> 
#> , , 2, 2
#> 
#>       I2
#> I1     [,1] [,2]
#>   [1,]    1    0
#>   [2,]    0    3
#> 
#> attr(,"class")
#> [1] "tensor"
```
