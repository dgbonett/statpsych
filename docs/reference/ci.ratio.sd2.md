# Confidence interval for a 2-group ratio of standard deviations

Computes a robust confidence interval for a ratio of population standard
deviations in a 2-group design. This function is a modification of the
confidence interval proposed by Bonett (2006). The original Bonett
method used a pooled kurtosis estimate in the standard error that
assumed equal variances, which limited the confidence interval's use to
tests of equal population variances and equivalence tests. This function
uses a pooled kurtosis estimate that does not assume equal variances and
provides a useful confidence interval for a ratio of standard deviations
under general conditions. This function requires of minimum sample size
of four per group but sample sizes of at least 10 per group are
recommended.

## Usage

``` r
ci.ratio.sd2(alpha, y1, y2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y1:

  vector of scores for group 1

- y2:

  vector of scores for group 2

## Value

Returns a 1-row matrix. The columns are:

- SD1 - estimated SD for group 1

- SD2 - estimated SD for group 2

- SD1/SD2 - estimate of SD ratio

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2006). “Robust confidence interval for a ratio of standard
deviations.” *Applied Psychological Measurement*, **30**(5), 432–439.
ISSN 0146-6216,
[doi:10.1177/0146621605279551](https://doi.org/10.1177/0146621605279551)
.

## Examples

``` r
y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
ci.ratio.sd2(.05, y1, y2)
#>       SD1      SD2   SD1/SD2       LL       UL
#>  5.711587 6.450667 0.8854257 0.486279 1.728396

# Should return:
#      SD1      SD2    SD1/SD2       LL       UL
# 5.711587 6.450667  0.8854257 0.486279 1.728396

```
