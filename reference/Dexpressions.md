# List of unevaluated expressions for (mixed) partial derivatives of fitness with respect to linear predictors.

Unevaluated expressions for (mixed) partial derivatives of fitness with
respect to linear predictors for survival and fecundity.

## Usage

``` r
Dexpressions
```

## Value

- PW.d0W:

  Fitness (W) function for the Poisson-Weibull (PW) model.

- PW.d1Wds:

  First Partial derivative of fitness (d1W) with respect to survival
  (d1s) linear predictor for the Poisson-Weibull (PW) model.

- PW.d1Wdf:

  First Partial derivative of fitness (d1W) with respect to fecundity
  (d1f) linear redictor for the Poisson-Weibull (PW) model.

- PW.d3Wd2sd1f:

  Mixed third partial derivative of fitness (d3W) with 2nd derivative of
  survival linear predictor (d2s) and first derivative of fecundity
  linear predictor (d1f) from the Poisson-Weibull (PW) model.

- PW.d3Wdsd2f:

  and so on ...

- PW.d2Wd2f:

- PW.d2Wd2s:

- PW.d3Wd3s:

- PW.d3Wd3f:

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`Dtensor`](https://jarrodhadfield.github.io/MCMCglmm/reference/Dtensor.md)
