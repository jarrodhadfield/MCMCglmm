# Estimates the marginal parameter modes using kernel density estimation

Estimates the marginal parameter modes using kernel density estimation

## Usage

``` r
posterior.mode(x, adjust=0.1, ...)
```

## Arguments

- x:

  mcmc object

- adjust:

  numeric, passed to [`density`](https://rdrr.io/r/stats/density.html)
  to adjust the bandwidth of the kernal density

- ...:

  other arguments to be passed

## Value

modes of the kernel density estimates

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## See also

[`density`](https://rdrr.io/r/stats/density.html)

## Examples

``` r
v<-rIW(as.matrix(1),10, n=1000)
hist(v)
abline(v=posterior.mode(mcmc(v)), col="red")
```
