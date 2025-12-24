# Summarising GLMM Fits from MCMCglmm

`summary` method for class `"MCMCglmm"`. The returned object is suitable
for printing with the `print.summary.MCMCglmm` method.

## Usage

``` r
# S3 method for class 'MCMCglmm'
summary(object, random=FALSE, ...)
```

## Arguments

- object:

  an object of class `"MCMCglmm"`

- random:

  logical: should the random effects be summarised

- ...:

  Further arguments to be passed

## Value

- DIC:

  Deviance Information Criterion

- fixed.formula:

  model formula for the fixed terms

- random.formula:

  model formula for the random terms

- residual.formula:

  model formula for the residual terms

- solutions:

  posterior mean, 95% HPD interval, MCMC p-values and effective sample
  size of fixed (and random) effects

- Gcovariances:

  posterior mean, 95% HPD interval and effective sample size of random
  effect (co)variance components

- Gterms:

  indexes random effect (co)variances by the component terms defined in
  the random formula

- Rcovariances:

  posterior mean, 95% HPD interval and effective sample size of residual
  (co)variance components

- Rterms:

  indexes residuals (co)variances by the component terms defined in the
  rcov formula

- csats:

  chain length, burn-in and thinning interval

- cutpoints:

  posterior mean, 95% HPD interval and effective sample size of
  cut-points from an ordinal model

- theta_scale:

  posterior mean, 95% HPD interval, MCMC p-values and effective sample
  size of scaling parameter in theta_scale models.

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`MCMCglmm`](https://jarrodhadfield.github.io/MCMCglmm/reference/MCMCglmm.md)
