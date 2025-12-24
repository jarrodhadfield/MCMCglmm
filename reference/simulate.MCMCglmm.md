# Simulate method for GLMMs fitted with MCMCglmm

Simulated response vectors for GLMMs fitted with MCMCglmm

## Usage

``` r
# S3 method for class 'MCMCglmm'
simulate(object, nsim = 1, seed = NULL, newdata=NULL, marginal = object$Random$formula, 
          type = "response", it=NULL, posterior = "all", verbose=FALSE, ...)
```

## Arguments

- object:

  an object of class `"MCMCglmm"`

- nsim:

  number of response vectors to simulate. Defaults to `1`.

- seed:

  Either `NULL` or an integer that will be used in a call to `set.seed`
  before simulating the response vectors. The default, `NULL` will not
  change the random generator state.

- newdata:

  An optional data frame for which to simulate new observations

- marginal:

  formula defining random effects to be maginalised

- type:

  character; either "terms" (link scale) or "response" (data scale)

- it:

  integer; optional, MCMC iteration on which predictions should be based

- posterior:

  character; if `it` is `NULL` should the response vector be simulated
  using the marginal posterior means ("mean") of the parameters, or the
  posterior modes ("mode"), random draws from the posterior with
  replacement ("distribution") or without replacement ("all")

- verbose:

  logical; if `TRUE`, warnings are issued with newdata when the original
  model has fixed effects that do not appear in newdata and/or newdata
  has random effects not present in the original model.

- ...:

  Further arguments to be passed

## Value

A matrix (with nsim columns) of simulated response vectors

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`MCMCglmm`](https://jarrodhadfield.github.io/MCMCglmm/reference/MCMCglmm.md)
