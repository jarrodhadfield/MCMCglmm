# Resolves a Prior Specification in MCMCglmm

Generates the (co)variance matrix prior parameters in MCMCglmm given a
prior specification that might utilise a `prior_generator`.

## Usage

``` r
resolve_prior(prior_element, k=NULL)
```

## Arguments

- prior_element:

  prior spefictaion for a (co)variance matrix

- k:

  dimension of the (co)variance matrix

- vtype:

  type of covariance structure used: currently only `us` and `idh` are
  allowed

## Details

Primarily used internally to generate the (co)variance matrix prior
parameters used in MCMCglmm: `V`, `nu`, `alpha.mu` and `alpha.V`. If
these are already passed as a list, the list is simply returned. If a
`prior_generator` is used then the prior parameters are generated given
dimension of the (co)variance matrix, `k` and its stucture type `vtype`.

## Value

list of prior parameters for MCMCglmm

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

`prior_generator`,[`IW`](https://jarrodhadfield.github.io/MCMCglmm/reference/prior_generators.md),
[`IG`](https://jarrodhadfield.github.io/MCMCglmm/reference/prior_generators.md),
[`F`](https://jarrodhadfield.github.io/MCMCglmm/reference/prior_generators.md),
[`tSD`](https://jarrodhadfield.github.io/MCMCglmm/reference/prior_generators.md)

## Examples

``` r
prior_element1<-list(V=diag(2)/2, nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*1000)

resolve_prior(prior_element1)
#> $V
#>      [,1] [,2]
#> [1,]  0.5  0.0
#> [2,]  0.0  0.5
#> 
#> $nu
#> [1] 2
#> 
#> $alpha.mu
#> [1] 0 0
#> 
#> $alpha.V
#>      [,1] [,2]
#> [1,] 1000    0
#> [2,]    0 1000
#> 

prior_element2<-F(df2=1, scale=1000)

resolve_prior(prior_element2, k=2, vtype="us")
#> $V
#>      [,1] [,2]
#> [1,]  0.5  0.0
#> [2,]  0.0  0.5
#> 
#> $nu
#> [1] 2
#> 
#> $alpha.mu
#> [1] 0 0
#> 
#> $alpha.V
#>      [,1] [,2]
#> [1,] 1000    0
#> [2,]    0 1000
#> 

resolve_prior(prior_element2, k=2, vtype="idh")
#> $V
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
#> $nu
#> [1] 1
#> 
#> $alpha.mu
#> [1] 0 0
#> 
#> $alpha.V
#>      [,1] [,2]
#> [1,] 1000    0
#> [2,]    0 1000
#> 
```
