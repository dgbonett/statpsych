# Confidence interval for Cramer's V

Computes a confidence interval for a population Cramer's V coefficient
of nominal association for an r x s contingency table. The confidence
interval is based on a noncentral chi-square distribution, and an
approximate standard error is recovered from the confidence interval.

For more details, see Section 3.10 of Bonett (2021, Volume 3)

## Usage

``` r
ci.cramer(alpha, chisqr, r, c, n)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- chisqr:

  Pearson chi-square test statistic of independence

- r:

  number of rows in contingency table

- c:

  number of columns in contingency table

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of Cramer's V

- SE - recovered standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Smithson M (2003). *Confidence Intervals*. Sage. Bonett DG (2021url).
*Statistical Methods for Psychologists*.

## Examples

``` r
ci.cramer(.05, 19.21, 2, 3, 200)
#>  Estimate     SE   LL    UL
#>      0.31 0.0718 0.16 0.442

# Should return:
# Estimate     SE   LL    UL
#     0.31 0.0718 0.16 0.442
 
```
