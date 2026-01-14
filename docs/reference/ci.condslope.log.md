# Confidence intervals for conditional (simple) slopes in a logistic model

Computes confidence intervals and test statistics for population
conditional slopes (simple slopes) in a logistic model that includes a
predictor variable (x1), a moderator variable (x2), and a product
predictor variable (x1\*x2). Conditional slopes are computed at low and
high values of the moderator variable.

For more details, see Section 4.9 of Bonett (2021, Volume 3)

## Usage

``` r
ci.condslope.log(alpha, b1, b2, se1, se2, cov, lo, hi)
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

## Value

Returns a 2-row matrix. The columns are:

- Estimate - estimated conditional slope

- exp(Estimate) - estimated exponentiated conditional slope

- z - z test statistic

- p - two-sided p-value

- LL - lower limit of the exponentiated confidence interval

- UL - upper limit of the exponentiated confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.condslope.log(.05, .132, .154, .031, .021, .015, 5.2, 10.6)
#>                   Estimate exp(Estimate)      z       p       LL        UL
#> At low moderator    0.9328      2.541616 2.2698 0.02322 1.135802  5.687444
#> At high moderator   1.7644      5.838068 2.9065 0.00365 1.776421 19.186357

# Should return:
#                   Estimate exp(Estimate)      z       p 
# At low moderator    0.9328      2.541616 2.2698 0.02322 
# At high moderator   1.7644      5.838068 2.9065 0.00365 
#                          LL        UL
# At low moderator   1.135802  5.687444
# At high moderator  1.776421 19.186357

```
