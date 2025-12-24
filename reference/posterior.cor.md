# Transforms posterior distribution of covariances into correlations

Transforms posterior distribution of covariances into correlations

## Usage

``` r
posterior.cor(x)
```

## Arguments

- x:

  mcmc object of (co)variances stacked column-wise

## Value

posterior correlation matrices

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`posterior.evals`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.evals.md),
[`posterior.inverse`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.inverse.md),
[`posterior.ante`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.ante.md)

## Examples

``` r
v<-rIW(diag(2),3, n=1000)
hist(posterior.cor(mcmc(v))[,2])
```
