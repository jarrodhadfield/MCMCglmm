# d-divergence

Calculates Ovaskainen's (2008) d-divergence between 2 zero-mean
multivariate normal distributions.

## Usage

``` r
Ddivergence(CA=NULL, CB=NULL, n=10000)
```

## Arguments

- CA:

  Matrix A

- CB:

  Matrix B

- n:

  number of Monte Carlo samples for approximating the integral

## Value

d-divergence

## Note

In versions of MCMCglmm \<2.26 Ovaskainen's (2008) d-divergence was
incorrectly calculated.

## References

Ovaskainen, O. et. al. (2008) Proc. Roy. Soc - B (275) 1635 593-750

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
CA<-rIW(diag(2),10, n=1)
CB<-rIW(diag(2),10, n=1)
Ddivergence(CA, CB)
#> [1] 0.2971437
```
