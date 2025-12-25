# Random Generation from the Conditional Inverse Wishart Distribution

Samples from the inverse Wishart distribution, with the possibility of
conditioning on a diagonal submatrix

## Usage

``` r
rIW(V, nu, fix=NULL, n=1, CM=NULL)
```

## Arguments

- V:

  Expected (co)varaince matrix as `nu` tends to infinity

- nu:

  degrees of freedom

- fix:

  optional integer indexing the partition to be conditioned on

- n:

  integer: number of samples to be drawn

- CM:

  matrix: optional matrix to condition on. If not given, and
  `fix!=NULL`, V_22 is conditioned on

## Details

If \\{\bf W^{-1}}\\ is a draw from the inverse Wishart, `fix` indexes
the diagonal element of \\{\bf W^{-1}}\\ which partitions \\{\bf
W^{-1}}\\ into 4 submatrices. `fix` indexes the upper left corner of the
lower diagonal matrix and it is this matrix that is conditioned on.

For example partioning \\{\bf W^{-1}}\\ such that

\$\$ {\bf W^{-1}} = \left\[ \begin{array}{cc} {\bf W^{-1}}\_{11}&{\bf
W^{-1}}\_{12}\\ {\bf W^{-1}}\_{21}&{\bf W^{-1}}\_{22}\\ \end{array}
\right\] \$\$ \$\$\$\$

fix indexes the upper left corner of \\{\bf W^{-1}}\_{22}\\. If
`CM!=NULL` then \\{\bf W^{-1}}\_{22}\\ is fixed at `CM`, otherwise
\\{\bf W^{-1}}\_{22}\\ is fixed at \\\texttt{V}\_{22}\\. For example, if
`dim(V)`=4 and `fix=2` then \\{\bf W^{-1}}\_{11}\\ is a 1X1 matrix and
\\{\bf W^{-1}}\_{22}\\ is a 3X3 matrix.

## Value

if `n` = 1 a matrix equal in dimension to `V`, if `n`\>1 a matrix of
dimension `n` x `length(V)`

## Note

In versions of MCMCglmm \>1.10 the arguments to `rIW` have changed so
that they are more intuitive in the context of
[`MCMCglmm`](https://jarrodhadfield.github.io/MCMCglmm/reference/MCMCglmm.md).
Following the notation of Wikipedia
(<https://en.wikipedia.org/wiki/Inverse-Wishart_distribution>) the
inverse scale matrix \\{\bm \Psi}=(\texttt{V\*nu})\\. In earlier
versions of MCMCglmm (\<1.11) \\{\bm \Psi} = \texttt{V}^{-1}\\. Although
the old parameterisation is consistent with the
[`riwish`](https://rdrr.io/pkg/MCMCpack/man/InvWishart.html) function in
MCMCpack and the
[`rwishart`](https://rdrr.io/pkg/bayesm/man/rwishart.html) function in
bayesm it is inconsistent with the prior definition for
[`MCMCglmm`](https://jarrodhadfield.github.io/MCMCglmm/reference/MCMCglmm.md).
The following pieces of code are sampling from the same distributions:

|                                                                                  |                       |
|----------------------------------------------------------------------------------|-----------------------|
| [`riwish`](https://rdrr.io/pkg/MCMCpack/man/InvWishart.html)`(nu, nu*V)`         | from MCMCpack         |
| [`rwishart`](https://rdrr.io/pkg/bayesm/man/rwishart.html)`(nu, solve(nu*V))$IW` | from bayesm           |
| `rIW(nu, solve(nu*V))`                                                           | from MCMCglmm \<1.11  |
| `rIW(V, nu)`                                                                     | from MCMCglmm \>=1.11 |

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## References

Korsgaard, I.R. et. al. 1999 Genetics Selection Evolution 31 (2) 177:181

## See also

[`rwishart`](https://rdrr.io/pkg/bayesm/man/rwishart.html),
[`rwish`](https://rdrr.io/pkg/MCMCpack/man/Wishart.html)

## Examples

``` r
nu<-10
V<-diag(4)
rIW(V, nu, fix=2)
#>             [,1]        [,2]      [,3]      [,4]
#> [1,]  0.71251554 -0.02395105 0.2857713 0.2408761
#> [2,] -0.02395105  1.00000000 0.0000000 0.0000000
#> [3,]  0.28577129  0.00000000 1.0000000 0.0000000
#> [4,]  0.24087612  0.00000000 0.0000000 1.0000000
```
