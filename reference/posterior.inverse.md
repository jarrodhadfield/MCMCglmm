# Posterior distribution of matrix inverse

Posterior distribution of matrix inverse

## Usage

``` r
posterior.inverse(x)
```

## Arguments

- x:

  mcmc object of (co)variances stacked column-wise

## Value

posterior of inverse (co)variance matrices

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`posterior.cor`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.cor.md),
[`posterior.evals`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.evals.md),
[`posterior.ante`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.ante.md)

## Examples

``` r
  V<-rIW(diag(2),3, n=1000)
  inv.V<-posterior.inverse(V)
#> Warning: posterior.inverse expecting mcmc object
```
