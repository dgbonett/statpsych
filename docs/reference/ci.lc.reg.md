# Confidence interval for a linear contrast of regression coefficients in multiple group regression model

Computes a confidence interval and test statistic for a linear contrast
of population regression coefficients (e.g., a y-intercept or a slope
coefficient) across groups in a multiple group regression model.
Equality of error variances across groups is not assumed. A
Satterthwaite adjustment to the degrees of freedom is used to improve
the accuracy of the confidence interval.

For more details, see Section 2.20 of Bonett (2021, Volume 2)

## Usage

``` r
ci.lc.reg(alpha, est, se, n, s, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- est:

  vector of parameter estimates

- se:

  vector of standard errors

- n:

  vector of group sample sizes

- s:

  number of predictor variables for each within-group model

- v:

  vector of contrast coefficients

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated linear contrast

- SE - standard error

- t - t test statistic

- df - degrees of freedom

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
est <- c(1.74, 1.83, 0.482)
se <- c(.483, .421, .395)
n <- c(40, 40, 40)
v <- c(.5, .5, -1)
ci.lc.reg(.05, est, se, n, 4, v)
#>  Estimate        SE     t      df       p        LL       UL
#>     1.303 0.5085838 2.562 78.8197 0.01231 0.2906532 2.315347

# Should return:
# Estimate        SE      t      df       p        LL       UL
#    1.303 0.5085838 2.5620 78.8197 0.01231 0.2906532 2.315347
 
```
