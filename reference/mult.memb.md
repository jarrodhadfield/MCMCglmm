# Design Matrices for Multiple Membership Models

Forms design matrices for multiple membership models

## Usage

``` r
mult.memb(formula)
```

## Arguments

- formula:

  formula

## Value

design matrix

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Details

Currently `mult.memb` can only usefully be used inside an `idv` variance
function. The formula usually contains serveral factors that have the
same factor levels.

## Examples

``` r
fac1<-factor(sample(letters[1:3], 5, TRUE), levels=letters[1:3])
fac2<-factor(sample(letters[1:3], 5, TRUE), levels=letters[1:3])
cbind(fac1, fac2)
#>      fac1 fac2
#> [1,]    2    3
#> [2,]    3    2
#> [3,]    2    2
#> [4,]    3    1
#> [5,]    1    2
mult.memb(~fac1+fac2)
#>   fac1a fac1b fac1c
#> 1     0     1     1
#> 2     0     1     1
#> 3     0     2     0
#> 4     1     0     1
#> 5     1     1     0
#> attr(,"assign")
#> [1] 1 1 1
#> attr(,"contrasts")
#> attr(,"contrasts")$fac1
#> [1] "contr.treatment"
#> 
```
