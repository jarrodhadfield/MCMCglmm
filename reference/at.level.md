# Incidence Matrix of Levels within a Factor

Incidence matrix of levels within a factor

## Usage

``` r
at.level(x, level)
```

## Arguments

- x:

  factor

- level:

  factor level

## Value

incidence matrix for level in x

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`at.set`](https://jarrodhadfield.github.io/MCMCglmm/reference/at.set.md)

## Examples

``` r
fac<-gl(3,10,30, labels=letters[1:3])
x<-rnorm(30)
model.matrix(~at.level(fac,"b"):x)
#>    (Intercept) at.level(fac, "b"):x
#> 1            1            0.0000000
#> 2            1            0.0000000
#> 3            1            0.0000000
#> 4            1            0.0000000
#> 5            1            0.0000000
#> 6            1            0.0000000
#> 7            1            0.0000000
#> 8            1            0.0000000
#> 9            1            0.0000000
#> 10           1            0.0000000
#> 11           1           -0.9561440
#> 12           1            2.0181869
#> 13           1           -1.3096437
#> 14           1           -0.1127342
#> 15           1           -0.2806937
#> 16           1           -0.5345137
#> 17           1            0.5748624
#> 18           1            1.2344475
#> 19           1            0.1281139
#> 20           1           -0.5060917
#> 21           1            0.0000000
#> 22           1            0.0000000
#> 23           1            0.0000000
#> 24           1            0.0000000
#> 25           1            0.0000000
#> 26           1            0.0000000
#> 27           1            0.0000000
#> 28           1            0.0000000
#> 29           1            0.0000000
#> 30           1            0.0000000
#> attr(,"assign")
#> [1] 0 1
```
