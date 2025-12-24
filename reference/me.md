# Design Matrix for Measurement Error Model

Sets up design matrix for measurement error models.

## Usage

``` r
me(formula, error=NULL, group=NULL, type="classical")
```

## Arguments

- formula:

  [`formula`](https://rdrr.io/r/stats/formula.html) for the fixed
  effects.

- error:

  character; name of column in
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) in which
  standard error (`type="classical"` or `type="berkson"`) or
  miscalssification error (`type="dclassical"`) is stored.

- group:

  name of column in
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) in which groups
  are stored. Rows of the design matrix with the same `group` level are
  assumed to pertain to the same obsevation of the covariate that is
  measured with error.

- type:

  character; one of `type="classical"`, `type="berkson"`,
  `type="dclassical"` or `type="dberkson"` (see details)

## Value

design matrix, with a prior distribution attribute

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>
