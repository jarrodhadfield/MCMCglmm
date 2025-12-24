# Tensor of (mixed) partial derivatives

Forms tensor of (mixed) partial derivatives

## Usage

``` r
Dtensor(expr, name=NULL, mu = NULL, m=1, evaluate = TRUE)
```

## Arguments

- expr:

  'expression'

- name:

  character vector, giving the variable names with respect to which
  derivatives will be computed. If NULL all variables in the expression
  will be used

- mu:

  optional: numeric vector, at which the derivatives are evaluated

- m:

  order of derivative

- evaluate:

  logical; if `TRUE` the derivatives are evaluated at `mu`, if `FALSE`
  the derivatives are left unevaluated

## Value

- Dtensor:

  (list) of unevaluated expression(s) if `evaluate=FALSE` or a tensor if
  `evaluate=TRUE`

## References

Rice, S.H. (2004) Evolutionary Theory: Mathematical and Conceptual
Foundations. Sinauer (MA) USA.

## Author

Jarrod Hadfield j.hadfield@ed.ac.uk

## See also

[`evalDtensor`](https://jarrodhadfield.github.io/MCMCglmm/reference/evalDtensor.md),
[`Dexpressions`](https://jarrodhadfield.github.io/MCMCglmm/reference/Dexpressions.md),
[`D`](https://rdrr.io/r/stats/deriv.html)

## Examples

``` r
f<-expression(beta_1 + time * beta_2 + u)
Dtensor(f,eval=FALSE)
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> beta_2
#> 
#> [[3]]
#> time
#> 
#> [[4]]
#> [1] 1
#> 
```
