# Random Generation from a Truncated Normal Distribution

Samples from the Truncated Normal Distribution

## Usage

``` r
rtnorm(n = 1, mean = 0, sd = 1, lower = -Inf, upper = Inf)
```

## Arguments

- n:

  integer: number of samples to be drawn

- mean:

  vector of means

- sd:

  vector of standard deviations

- lower:

  left truncation point

- upper:

  right truncation point

## Value

vector

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## References

Robert, C.P. (1995) Statistics & Computing 5 121-125

## See also

[`rtnorm`](https://chjackson.github.io/msm/reference/tnorm.html)

## Examples

``` r
hist(rtnorm(100, lower=-1, upper=1))
```
