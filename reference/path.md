# Design Matrix for Path Analyses

Forms design matrix for path analyses that involve paths within residual
blocks

## Usage

``` r
path(cause=NULL, effect=NULL, k)
```

## Arguments

- cause:

  integer; index of predictor \`trait' within residual block

- effect:

  integer; index of response \`trait' within residual block

- k:

  integer; dimension of residual block

## Value

design matrix

## Note

For more general path anlaytic models see
[sir](https://jarrodhadfield.github.io/MCMCglmm/reference/sir.md) which
allows paths to exist between responses that are not in the same
residual block. However, `sir` does not handle non-Gaussian or missing
responses. Note that path models involving non-Gaussian data are defined
on the link scale which may not always be appropriate.

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[sir](https://jarrodhadfield.github.io/MCMCglmm/reference/sir.md)

## Examples

``` r
path(1, 2,2)
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    1    0
```
