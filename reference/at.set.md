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
#> 11           1         0.45812074
#> 12           1         0.28437239
#> 13           1         0.15192926
#> 14           1        -0.78147093
#> 15           1        -0.07033439
#> 16           1         0.43091323
#> 17           1        -0.76990282
#> 18           1         1.77283669
#> 19           1        -0.86333630
#> 20           1         0.41221574
#> 21           1         1.79032792
#> 22           1         0.63244744
#> 23           1        -0.47924247
#> 24           1        -0.33278895
#> 25           1         1.74818948
#> 26           1         0.92827136
#> 27           1         0.85546100
#> 28           1        -1.36348797
#> 29           1         0.07556184
#> 30           1        -1.71581037
#> attr(,"assign")
#> [1] 0 1
```
