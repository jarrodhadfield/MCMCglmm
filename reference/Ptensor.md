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
#>   [1,] 1.04358324 0.02483401
#>   [2,] 0.02483401 1.09542565
#> attr(,"class")
#> [1] "tensor" "matrix"
cov(y)*((n-1)/n)
#>            [,1]       [,2]
#> [1,] 1.04462891 0.02485889
#> [2,] 0.02485889 1.09652327
```
