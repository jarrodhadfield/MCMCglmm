# Forms expected (co)variances for GLMMs fitted with MCMCglmm

Forms the expected covariance structure of link-scale observations for
GLMMs fitted with MCMCglmm

## Usage

``` r
buildV(object, marginal=object$Random$formula, diag=TRUE, it=NULL, posterior="mean", ...)
```

## Arguments

- object:

  an object of class `"MCMCglmm"`

- marginal:

  formula defining random effects to be maginalised

- diag:

  logical; if `TRUE` the covariances betwween observations are not
  calculated

- it:

  integer; optional, MCMC iteration on which covariance structure should
  be based

- posterior:

  character; if `it` is `NULL` should the covariance structure be based
  on the marginal posterior means (`'mean'`) of the VCV parameters, or
  the posterior modes (`'mode'`), or a random draw from the posterior
  with replacement (`'distribution'`). If `posterior=="all"` the
  posterior distribution of observation variances is returned

- ...:

  Further arguments to be passed

## Value

If `diag=TRUE` an n by n covariance matrix. If `diag=FALSE` and
`posterior!="all"` an 1 by n matrix of variances. If `posterior=="all"`
an nit by n matrix of variances (where nit is the number of saved MCMC
iterations).

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`MCMCglmm`](https://jarrodhadfield.github.io/MCMCglmm/reference/MCMCglmm.md)
