# Posterior distribution of ante-dependence parameters

Posterior distribution of ante-dependence parameters

## Usage

``` r
posterior.ante(x,k=1)
```

## Arguments

- x:

  mcmc object of (co)variances stacked column-wise

- k:

  order of the ante-dependence structure

## Value

posterior ante-dependence parameters (innovation variances followed by
regression ceofficients)

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`posterior.cor`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.cor.md),
[`posterior.evals`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.evals.md),
[`posterior.inverse`](https://jarrodhadfield.github.io/MCMCglmm/reference/posterior.inverse.md)

## Examples

``` r
v<-rIW(diag(2),10, n=1000)
plot(posterior.ante(mcmc(v),1))
#> Error in colnames(ante)[1:k] <- colnames(x)[seq(1, k^2, k + 1)]: replacement has length zero
```
