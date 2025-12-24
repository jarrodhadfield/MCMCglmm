# Posterior distribution of eigenvalues

Posterior distribution of eigenvalues

## Usage

``` r
posterior.evals(x)
```

## Arguments

- x:

  mcmc object of (co)variances stacked column-wise

## Value

posterior eigenvalues

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`posterior.cor`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.cor.md),
[`posterior.inverse`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.inverse.md),
[`posterior.ante`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.ante.md)

## Examples

``` r
v<-rIW(diag(2),3, n=1000)
hist(posterior.evals(mcmc(v))[,2])
```
