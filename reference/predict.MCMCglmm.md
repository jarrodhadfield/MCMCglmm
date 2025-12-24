# Predict method for GLMMs fitted with MCMCglmm

Predicted values for GLMMs fitted with MCMCglmm

## Usage

``` r
# S3 method for class 'MCMCglmm'
predict(object, newdata=NULL, marginal=object$Random$formula,
        type="response", interval="none", level=0.95, it=NULL, 
        posterior="all", verbose=FALSE, approx="numerical", ...)
```

## Arguments

- object:

  an object of class `"MCMCglmm"`

- newdata:

  An optional data frame in which to look for variables with which to
  predict

- marginal:

  formula defining random effects to be maginalised

- type:

  character; either "terms" (link scale) or "response" (data scale)

- interval:

  character; either "none", "confidence" or "prediction"

- level:

  A numeric scalar in the interval (0,1) giving the target probability
  content of the intervals.

- it:

  integer; optional, MCMC iteration on which predictions should be based

- posterior:

  character; should marginal posterior predictions be calculated
  ("all"), or should they be made conditional on the marginal posterior
  means ("mean") of the parameters, the posterior modes ("mode"), or a
  random draw from the posterior ("distribution").

- verbose:

  logical; if `TRUE`, warnings are issued with newdata when the original
  model has fixed effects that do not appear in newdata and/or newdata
  has random effects not present in the original model.

- approx:

  character; for distributions for which the mean cannot be calculated
  analytically what approximation should be used: numerical integration
  (`numerical`; slow), second order Taylor expansion (`taylor2`) and for
  logistic models approximations presented in Diggle (2004) (`diggle`)
  and McCulloch and Searle (2001) (`mcculloch`)

- ...:

  Further arguments to be passed

## Value

Expectation and credible interval

## References

Diggle P, et al. (2004). Analysis of Longitudinal Data. 2nd Edition.
Oxford University Press.

McCulloch CE and Searle SR (2001). Generalized, Linear and Mixed Models.
John Wiley & Sons, New York.

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`MCMCglmm`](https://jarrodhadfield.github.io/MCMCglmm/reference/MCMCglmm.md)
