# Confidence interval for a linear contrast of means in a between-subjects design

Computes a test statistic and confidence interval for a linear contrast
of means in a between-subjects design. This function computes both
unequal variance and equal variance confidence intervals and test
statistics. A Satterthwaite adjustment to the degrees of freedom is used
with the unequal variance method.

For more details, see Section 3.3 of Bonett (2021, Volume 1)

## Usage

``` r
ci.lc.mean.bs(alpha, m, sd, n, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m:

  vector of estimated group means

- sd:

  vector of estimated group standard deviations

- n:

  vector of sample sizes

- v:

  vector of between-subjects contrast coefficients

## Value

Returns a 2-row matrix. The columns are:

- Estimate - estimated linear contrast

- SE - standard error

- t - t test statistic

- df - degrees of freedom

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
m <- c(24.9, 23.1, 16.4)
sd <- c(5.21, 4.67, 4.98)
n <- c(30, 30, 30)
v <- c(.5, .5, -1)
ci.lc.mean.bs(.05, m, sd, n, v)
#>                              Estimate       SE      t    df p       LL       UL
#> Equal Variances Assumed:          7.6 1.108703 6.8549 87.00 0 5.396332 9.803668
#> Equal Variances Not Assumed:      7.6 1.111135 6.8399 57.59 0 5.375484 9.824516

# Should return:
#                              Estimate       SE      t    df p       LL       UL
# Equal Variances Assumed:          7.6 1.108703 6.8549 87.00 0 5.396332 9.803668
# Equal Variances Not Assumed:      7.6 1.111135 6.8399 57.59 0 5.375484 9.824516

m <- c(33.5, 37.9, 38.0, 44.1)
sd <- c(3.84, 3.84, 3.65, 4.98)
n <- c(10, 10, 10, 10)
v <- c(.5, .5, -.5, -.5)
ci.lc.mean.bs(.05, m, sd, n, v)
#>                              Estimate       SE      t    df       p        LL
#> Equal Variances Assumed:        -5.35 1.300136 -4.115 36.00 0.00022 -7.986797
#> Equal Variances Not Assumed:    -5.35 1.300136 -4.115 33.52 0.00024 -7.993588
#>                                     UL
#> Equal Variances Assumed:     -2.713203
#> Equal Variances Not Assumed: -2.706412

# Should return:
#                              Estimate       SE      t    df       p 
# Equal Variances Assumed:        -5.35 1.300136 -4.115 36.00 0.00022 
# Equal Variances Not Assumed:    -5.35 1.300136 -4.115 33.52 0.00024 
#                                     LL        UL
# Equal Variances Assumed:     -7.986797 -2.713203
# Equal Variances Not Assumed: -7.993583 -2.706417

```
