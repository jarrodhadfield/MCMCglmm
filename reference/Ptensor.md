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
#>   [1,] 1.05107004 0.02747103
#>   [2,] 0.02747103 1.08917803
#> attr(,"class")
#> [1] "tensor" "matrix"
cov(y)*((n-1)/n)
#>            [,1]       [,2]
#> [1,] 1.05212322 0.02749855
#> [2,] 0.02749855 1.09026939
```
