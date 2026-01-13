# Tensor of Sample (Mixed) Central Moments

Forms a tensor of sample (mixed) central moments

## Usage

``` r
Ptensor(x, k)
```

## Arguments

- x:

  matrix; traits in columns samples in rows

- k:

  kth central moment

## Value

tensor

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
n<-1000
y<-matrix(rnorm(n), n/2, 2)
Ptensor(y,2)
#>       I2
#> I1           [,1]       [,2]
#>   [1,] 0.96135518 0.02785343
#>   [2,] 0.02785343 1.01490129
#> attr(,"class")
#> [1] "tensor" "matrix"
cov(y)*((n-1)/n)
#>            [,1]       [,2]
#> [1,] 0.96231846 0.02788134
#> [2,] 0.02788134 1.01591822
```
