# Incidence Matrix of Combined Levels within a Factor

Incidence Matrix of Combined Levels within a Factor

## Usage

``` r
at.set(x, level)
```

## Arguments

- x:

  factor

- level:

  set of factor levels

## Value

incidence matrix for the set `level` in x

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`at.level`](https://jarrodhadfield.github.io/MCMCglmm/reference/at.level.md)

## Examples

``` r
fac<-gl(3,10,30, labels=letters[1:3])
x<-rnorm(30)
model.matrix(~at.set(fac,2:3):x)
#>    (Intercept) at.set(fac, 2:3):x
#> 1            1          0.0000000
#> 2            1          0.0000000
#> 3            1          0.0000000
#> 4            1          0.0000000
#> 5            1          0.0000000
#> 6            1          0.0000000
#> 7            1          0.0000000
#> 8            1          0.0000000
#> 9            1          0.0000000
#> 10           1          0.0000000
#> 11           1         -0.6127966
#> 12           1         -0.2828786
#> 13           1          0.4892402
#> 14           1          0.1766764
#> 15           1          0.3064992
#> 16           1          1.4139759
#> 17           1         -0.4980426
#> 18           1         -0.5556104
#> 19           1          0.7263960
#> 20           1          0.1712834
#> 21           1         -1.1799395
#> 22           1          0.7144402
#> 23           1          0.5739545
#> 24           1          0.1679005
#> 25           1         -0.9858416
#> 26           1          1.3334297
#> 27           1         -0.6417486
#> 28           1          0.2688340
#> 29           1          2.0230553
#> 30           1          0.5236421
#> attr(,"assign")
#> [1] 0 1
```
