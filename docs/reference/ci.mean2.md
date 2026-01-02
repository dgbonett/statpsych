# Confidence interval for a 2-group mean difference

Computes equal variance and unequal variance confidence intervals for a
population 2-group mean difference using the estimated means, estimated
standard deviations, and sample sizes. Also computes equal variance and
unequal variance independent-samples t-tests. Use the t.test function
for raw data input.

For more details, see Section 2.3 of Bonett (2021, Volume 1)

## Usage

``` r
ci.mean2(alpha, m1, m2, sd1, sd2, n1, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  estimated mean for group 1

- m2:

  estimated mean for group 2

- sd1:

  estimated standard deviation for group 1

- sd2:

  estimated standard deviation for group 2

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a 2-row matrix. The columns are:

- Estimate - estimated mean difference

- SE - standard error

- t - t test statistic

- df - degrees of freedom

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa. Bonett DG (2021url). *Statistical Methods
for Psychologists*.

## Examples

``` r
ci.mean2(.05, 19.4, 11.3, 2.70, 2.10, 40, 40)
#>                              Estimate        SE       t    df p       LL
#> Equal Variances Assumed:          8.1 0.5408327 14.9769 78.00 0 7.023285
#> Equal Variances Not Assumed:      8.1 0.5408327 14.9769 73.54 0 7.022255
#>                                    UL
#> Equal Variances Assumed:     9.176715
#> Equal Variances Not Assumed: 9.177745
# Should return:
#                              Estimate        SE       t    df p
# Equal Variances Assumed:          8.1 0.5408327 14.9769 78.00 0
# Equal Variances Not Assumed:      8.1 0.5408327 14.9769 73.54 0
#                                    LL       UL
# Equal Variances Assumed:     7.023285 9.176715
# Equal Variances Not Assumed: 7.022256 9.177744

```
