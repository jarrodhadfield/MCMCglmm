# Generator Functions for Priors in MCMCglmm

Functions for generating (co)variance matrix prior specifications in
MCMCglmm that result in specified inverse-Wishart, inverse-gamma or
central-\\F\\ marginal priors for the variances.

## Usage

``` r
IW(V=1, nu=0.002)

  IG(shape=0.001, scale=0.001)
  
  F(df2=1, scale=1000)
  
  tSD(df=1, scale=sqrt(1000))
```

## Arguments

- V:

  expected varaince as `nu` tends to infinity in a scalar
  inverse-Wishart prior

- nu:

  degrees of freedom in a scalar inverse-Wishart prior

- shape:

  shape parameter of the inverse-gamma prior

- scale:

  scale parameter of the inverse-gamma, the scaled-F or scaled half-t

- df:

  degrees of freedom for the half-t prior on the standard deviation

- df2:

  denominator degrees of freedom for F prior (numerator
  degree-of-freedom is one)

## Details

Each genertor function returns a function that generates a list of prior
arguments need to specific (co)variance matrix priors in `MCMCglmm`.
Those prior arguments result in the marginal distributions for the
variancess being those specified in generator function. Since the
appropriate prior arguments depend on the dimension of the (co)variance
matrix, they are evalauted at run time once the dimension is determined
using `resove_prior`.

## Value

function of class `prior_generator`

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

`prior_generator`

## Examples

``` r
resolve_prior(F(df2=1, scale=1000), k=2)
#> Error in resolve_prior(F(df2 = 1, scale = 1000), k = 2): vtype must be specified
```
