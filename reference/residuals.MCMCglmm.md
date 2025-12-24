# Residuals form a GLMM fitted with MCMCglmm

`residuals` method for class `"MCMCglmm"`.

## Usage

``` r
# S3 method for class 'MCMCglmm'
residuals(object, type = c("deviance", "pearson", "working",
                                "response", "partial"), ...)
```

## Arguments

- object:

  an object of class `"MCMCglmm"`

- type:

  the type of residuals which should be returned. The alternatives are:
  `"deviance"` (default), `"pearson"`,`"working"`, `"response"`, and
  `"partial"`.

- ...:

  Further arguments to be passed

## Value

vector of residuals

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`residuals`](https://rdrr.io/r/stats/residuals.html),
[`MCMCglmm`](https://jarrodhadfield.github.io/MCMCglmm/reference/MCMCglmm.md)
