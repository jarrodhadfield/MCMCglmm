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
#> 1            1         0.00000000
#> 2            1         0.00000000
#> 3            1         0.00000000
#> 4            1         0.00000000
#> 5            1         0.00000000
#> 6            1         0.00000000
#> 7            1         0.00000000
#> 8            1         0.00000000
#> 9            1         0.00000000
#> 10           1         0.00000000
#> 11           1        -1.51357086
#> 12           1         0.83509548
#> 13           1         0.68736297
#> 14           1         0.39263047
#> 15           1         0.70734808
#> 16           1        -0.64231975
#> 17           1        -0.59173667
#> 18           1         0.03372468
#> 19           1        -0.84785405
#> 20           1        -2.34227823
#> 21           1        -0.43006191
#> 22           1         0.97764902
#> 23           1         0.07102241
#> 24           1        -0.35144927
#> 25           1         1.05181769
#> 26           1         0.80545320
#> 27           1        -0.12377926
#> 28           1        -0.74465437
#> 29           1         0.32145347
#> 30           1        -0.50822713
#> attr(,"assign")
#> [1] 0 1
```
