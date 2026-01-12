# Scheffe confidence interval for a linear contrast of means in a between-subjects design

Computes a Scheffe confidence interval for a linear contrast of
population means in a between-subjects design. A Scheffe p-value is
computed for the test statistic. The Scheffe method assumes equal
population variances. This function is useful in exploratory studies
where the linear contrast of means was not planned but was suggested by
the pattern of sample means. Use the
[ci.lc.mean.bs](https://dgbonett.github.io/statpsych/reference/ci.lc.mean.bs.md)
function with a Bonferroni adjusted alpha value to compute simultaneous
confidence intervals for two or more planned linear contrasts of means.

## Usage

``` r
ci.lc.mean.scheffe(alpha, m, sd, n, v)
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

- p - two-sided Scheffe p-value

- LL - lower limit of the Scheffe confidence interval

- UL - upper limit of the Scheffe confidence interval

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa.

## Examples

``` r
m <- c(33.5, 37.9, 38.0, 44.1)
sd <- c(3.49, 3.84, 3.65, 4.98)
n <- c(10, 10, 10, 10)
v <- c(.5, .5, -.5, -.5)
ci.lc.mean.scheffe(.05, m, sd, n, v)
#>  Estimate       SE       t       p        LL        UL
#>     -5.35 1.275231 -4.1953 0.00228 -9.089451 -1.610549

# Should return:
#  Estimate       SE       t       p        LL        UL
#     -5.35 1.275231 -4.1953 0.00228 -9.089451 -1.610549

```
