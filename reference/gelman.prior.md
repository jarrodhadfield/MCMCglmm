# Prior Covariance Matrix for Fixed Effects.

Prior covariance matrix for fixed effects had inputs been standardised
as suggested in Gelman et al. (2008).

## Usage

``` r
gelman.prior(formula, data, coef.scale=1, intercept.scale=coef.scale, singular.ok=FALSE)
```

## Arguments

- formula:

  [`formula`](https://rdrr.io/r/stats/formula.html) for the fixed
  effects.

- data:

  [`data.frame`](https://rdrr.io/r/base/data.frame.html).

- coef.scale:

  prior standard deviation for regression parameters (had inputs been
  standardised): default=1.

- intercept.scale:

  prior standard deviation for the intercept (had inputs been
  standardised): default=`coef.scale`.

- singular.ok:

  logical: if `FALSE` linear dependencies in the fixed effects are
  removed. if `TRUE` they are left in an estimated, although all
  information comes form the prior

- ...:

  Further arguments to be passed

## Value

list with elements `mu` (prior mean) and `V` prior covariance covariance
matrix

## References

Gelman, A. et al. (2008) The Annals of Appled Statistics 2 4 1360-1383

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Details

Gelman et al. (2008) suggest that the input variables in logistic
regression are standardised and that the associated regression
parameters are assumed independent in the prior. Gelman et al. (2008)
recommend a scaled t-distribution with a single degree of freedom
(scaled Cauchy) and a scale of 10 for the intercept and 2.5 for the
regression parameters. If the degree of freedom is infinity (i.e. a
normal distribution) then a prior covariance matrix `B$V` can be defined
for the regression parameters without input standardisation that
corresponds to a diagonal prior \\{\bf D}\\ for the regression
parameters had the inputs been standardised. The diagonal elements of
\\{\bf D}\\ are set to `coef.scale^2` except the first which is set to
`intercept.scale^2`. Depending on the link-function, the presence of
random effects, and the strength of prior required, suitable values for
`coef.scale` and `intercept.scale` may differ than those recommened by
Gelman et al. (2008). For details see
<https://jarrodhadfield.github.io/MCMCglmm/course-notes/glm.html#gelman-prior-sec>.
