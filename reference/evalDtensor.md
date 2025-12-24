# Evaluates a list of (mixed) partial derivatives

Evaluates a list of (mixed) partial derivatives

## Usage

``` r
evalDtensor(x, mu, m=1)
```

## Arguments

- x:

  unevaluated (list) of expression(s)

- mu:

  values at which the derivatives are evaluated: names need to match
  terms in x

- m:

  order of derivative

## Value

tensor

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`Dtensor`](https://jarrodhadfield.github.io/MCMCglmm/reference/Dtensor.md),
[`D`](https://rdrr.io/r/stats/deriv.html)

## Examples

``` r
f<-expression(beta_1 + time*beta_2+u)
Df<-Dtensor(f, eval=FALSE, m=2)
evalDtensor(Df, mu=data.frame(beta_1=0.5, beta_2=1, time=3, u=2.3))
#>  [1] 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0
#> attr(,"class")
#> [1] "tensor"
Dtensor(f, mu=c(1,3,1,2.3), m=2)
#>       I2
#> I1     [,1] [,2] [,3] [,4]
#>   [1,]    0    0    0    0
#>   [2,]    0    0    1    0
#>   [3,]    0    1    0    0
#>   [4,]    0    0    0    0
#> attr(,"class")
#> [1] "tensor" "matrix"
```
