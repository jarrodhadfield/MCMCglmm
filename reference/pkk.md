# Probability that all multinomial categories have a non-zero count.

Calculates the probability that all categories in a multinomial have a
non-zero count.

## Usage

``` r
pkk(prob, size)
```

## Arguments

- prob:

  numeric non-negative vector of length K, specifying the probability
  for the K classes; is internally normalized to sum 1. Infinite and
  missing values are not allowed.

- size:

  integer, say N, specifying the total number of objects that are put
  into K boxes in the typical multinomial experiment.

## Value

probability that there is at least one object in each of the K boxes

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
p<-runif(4)
pkk(p, 10)
#> [1] 0.1656116
```
