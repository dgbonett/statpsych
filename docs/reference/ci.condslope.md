# Confidence intervals for conditional (simple) slopes in a linear model

Computes confidence intervals and test statistics for population
conditional slopes (simple slopes) in a general linear model that
includes a predictor variable (x1), a moderator variable (x2), and a
product predictor variable (x1\*x2). Conditional slopes are computed at
specified low and high values of the moderator variable.

For more details, see Section 2.13 of Bonett (2021, Volume 2)

## Usage

``` r
ci.condslope(alpha, b1, b2, se1, se2, cov, lo, hi, dfe)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- b1:

  estimated slope coefficient for predictor variable

- b2:

  estimated slope coefficient for product variable

- se1:

  standard error for predictor coefficient

- se2:

  standard error for product coefficient

- cov:

  estimated covariance between predictor and product coefficients

- lo:

  low value of moderator variable

- hi:

  high value of moderator variable

- dfe:

  error degrees of freedom

## Value

Returns a 2-row matrix. The columns are:

- Estimate - estimated conditional slope

- t - t test statistic

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.condslope(.05, .132, .154, .031, .021, .015, 5.2, 10.6, 122)
#>                   Estimate        SE        t  df           p        LL
#> At low moderator    0.9328 0.4109570 2.269824 122 0.024973618 0.1192696
#> At high moderator   1.7644 0.6070517 2.906507 122 0.004342076 0.5626805
#>                         UL
#> At low moderator  1.746330
#> At high moderator 2.966119

# Should return:
#                   Estimate        SE        t  df           p 
# At low moderator    0.9328 0.4109570 2.269824 122 0.024973618 
# At high moderator   1.7644 0.6070517 2.906507 122 0.004342076 
#                           LL       UL
# At low moderator   0.1192696 1.746330
# At high moderator  0.5626805 2.966119
 
```
