# Design Matrix for Simultaneous and Recursive Relationships between Responses

Forms design matrix for simultaneous and recursive relationships between
responses

## Usage

``` r
sir(formula1=NULL, formula2=NULL, diag0=FALSE)
```

## Arguments

- formula1:

  formula

- formula2:

  formula

- diag0:

  logical: should the design matrix have zero's along the diagonal

## Value

design matrix

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
fac1<-factor(sample(letters[1:3], 5, TRUE), levels=letters[1:3])
fac2<-factor(sample(letters[1:3], 5, TRUE), levels=letters[1:3])
cbind(fac1, fac2)
#>      fac1 fac2
#> [1,]    3    3
#> [2,]    3    1
#> [3,]    1    2
#> [4,]    2    1
#> [5,]    3    2
sir(~fac1, ~fac2)
#>   1 2 3 4 5
#> 1 1 0 0 0 0
#> 2 1 0 0 0 0
#> 3 0 1 0 1 0
#> 4 0 0 1 0 1
#> 5 1 0 0 0 0
```
