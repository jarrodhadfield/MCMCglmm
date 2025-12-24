# Central Moments of a Uniform Distribution

Returns the central moments of a uniform distribution

## Usage

``` r
kunif(min, max, k)
```

## Arguments

- min, max:

  lower and upper limits of the distribution. Must be finite.

- k:

  k central moment, must be even

## Value

kth central moment

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`dunif`](https://rdrr.io/r/stats/Uniform.html)

## Examples

``` r
kunif(-1,1,4)
#> [1] 0.2
y<-runif(1000,-1,1)
mean((y-mean(y))^4)
#> [1] 0.1943839
```
