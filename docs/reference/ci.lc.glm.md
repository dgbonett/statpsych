# Confidence interval for a linear contrast of general linear model parameters

Computes the estimate, standard error, and confidence interval for a
linear contrast of parameters in a general linear model using
coef(object) and vcov(object) where "object" is a fitted model object
from the lm function.

For more details, see Section 2.28 of Bonett (2021, Volume 2)

## Usage

``` r
ci.lc.glm(alpha, n, b, V, q)
```

## Arguments

- alpha:

  alpha for 1 - alpha confidence

- n:

  sample size

- b:

  vector of parameter estimates from coef(object)

- V:

  covariance matrix of parameter estimates from vcov(object)

- q:

  vector of contrast coefficients

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of linear function

- SE - standard error

- t - t test statistic

- df - degrees of freedom

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
y <- c(43, 62, 49, 60, 36, 79, 55, 42, 67, 50)
x1 <- c(3, 6, 4, 6, 2, 7, 4, 2, 7, 5)
x2 <- c(4, 6, 3, 7, 1, 9, 3, 3, 8, 4)
out <- lm(y ~ x1 + x2)
b <- coef(out)
V <- vcov(out)
n <- length(y)
q <- c(0, .5, .5)
b
#> (Intercept)          x1          x2 
#>   26.891111    3.648889    2.213333 
ci.lc.glm(.05, n, b, V, q)
#>  Estimate        SE      t df       p       LL       UL
#>  2.931111 0.4462518 6.5683  7 0.00031 1.875893 3.986329

#  Should return:
# (Intercept)          x1          x2 
#   26.891111    3.648889    2.213333 
#
# Estimate        SE      t df       p       LL       UL
# 2.931111 0.4462518 6.5683  7 0.00031 1.875893 3.986329

```
