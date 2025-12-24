# Plots covariance matrices

Represents covariance matrices as 3-d ellipsoids using the `rgl`
package. Covariance matrices of dimension greater than 3 are plotted on
the subspace defined by the first three eigenvectors.

## Usage

``` r
plotsubspace(CA, CB=NULL, corr = FALSE, shadeCA = TRUE, 
      shadeCB = TRUE, axes.lab = FALSE, ...)
```

## Arguments

- CA:

  Matrix

- CB:

  Optional second matrix

- corr:

  If `TRUE` the covariance matrices are transformed into correlation
  matrices

- shadeCA:

  If `TRUE` the ellipsoid is solid, if `FALSE` the ellipsoid is
  wireframe

- shadeCB:

  If `TRUE` the ellipsoid is solid, if `FALSE` the ellipsoid is
  wireframe

- axes.lab:

  If `TRUE` the axes are labelled with the eigenvectors

- ...:

  further arguments to be passed

## Details

The matrix CA is always red, and the matrix CB if given is always blue.
The subspace is defined by the first three eigenvectors of CA, and the
percentage of variance for each matrix along these three dimensions is
given in the plot title.

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk> with code taken from the `rgl`
package

## See also

[`rgl`](https://dmurdoch.github.io/rgl/dev/reference/rgl-package.html)

## Examples

``` r
 if(requireNamespace("rgl")!=FALSE){
   G1<-rIW(diag(4),10)
   G2<-G1*1.2
 #  plotsubspace(G1, G2, shadeCB=FALSE)
 # commented out because of problems with rgl 
 } 
#> Loading required namespace: rgl
#> Warning: RGL: unable to open X11 display
#> Warning: 'rgl.init' failed, will use the null device.
#> See '?rgl.useNULL' for ways to avoid this warning.
```
