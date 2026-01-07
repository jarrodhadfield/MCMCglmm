# Prior Covariance Matrix for Fixed Effects.

Prior Covariance Matrix for Fixed Effects.

## Usage

``` r
gelman.prior(formula, data, scale=1, intercept=scale, singular.ok=FALSE)
```

## Arguments

- formula:

  [`formula`](https://rdrr.io/r/stats/formula.html) for the fixed
  effects.

- data:

  [`data.frame`](https://rdrr.io/r/base/data.frame.html).

- intercept:

  prior standard deviation for the intercept

- scale:

  prior standard deviation for regression parameters

- singular.ok:

  logical: if `FALSE` linear dependencies in the fixed effects are
  removed. if `TRUE` they are left in an estimated, although all
  information comes form the prior

## Value

prior covariance matrix

## References

Gelman, A. et al. (2008) The Annals of Appled Statistics 2 4 1360-1383

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Details

Gelman et al. (2008) suggest that the input variables of a categorical
regression are standardised and that the associated regression
parameters are assumed independent in the prior. Gelman et al. (2008)
recommend a scaled t-distribution with a single degree of freedom
(scaled Cauchy) and a scale of 10 for the intercept and 2.5 for the
regression parameters. If the degree of freedom is infinity (i.e. a
normal distribution) then a prior covariance matrix `B$V` can be defined
for the regression parameters without input standardisation that
corresponds to a diagonal prior \\{\bf D}\\ for the regression
parameters had the inputs been standardised. The diagonal elements of
\\{\bf D}\\ are set to `scale^2` except the first which is set to
`intercept^2`. With `family="binomial"`, \\D=\pi^{2}/3+\sigma^{2}\\
gives a prior that is approximately flat on the probability scale, where
\\\sigma^{2}\\ is the total variance due to the random effects and
residuals. With `family="threshold"` it is \\D=\sigma^{2}\\ and for
`family="ordinal"` (now superseded by `"threshold"`) it is
\\D=1+\sigma^{2}\\.

## Examples

``` r
dat<-data.frame(y=c(0,0,1,1), x=gl(2,2))
# data with complete separation

##############################
# standard probit regression #
##############################

prior1<-list(
  B=list(mu=c(0,0), V=gelman.prior(~x, data=dat, scale=1)), 
  R=list(V=1,fix=1))

m1<-MCMCglmm(y~x, prior=prior1, data=dat, family="threshold", verbose=FALSE)

p1<-pnorm(m1$Sol[,1]) # marginal probability when x=1

##################################################
# logistic regression with residual variance = 1 #
##################################################

prior2<-list(B=list(mu=c(0,0), V=gelman.prior(~x, data=dat, scale=sqrt(1+pi^2/3))),
             R=list(V=1,fix=1))

m2<-MCMCglmm(y~x, prior=prior2, data=dat, family="categorical", verbose=FALSE)

c2 <- (16 * sqrt(3)/(15 * pi))^2
p2<-plogis(m2$Sol[,1]/sqrt(1+c2)) # marginal probability when x=1

plot(mcmc.list(p1,p2))


```
