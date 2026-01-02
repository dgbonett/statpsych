# Test of a monotonic trend in medians for an ordered between-subjects factor

Computes simultaneous confidence intervals for all adjacent pairwise
comparisons of population medians using sample group medians and
standard errors as input. If one or more lower limits are greater than 0
and no upper limit is less than 0, then conclude that the population
medians are monotonic decreasing. If one or more upper limits are less
than 0 and no lower limits are greater than 0, then conclude that the
population medians are monotonic increasing. Reject the hypothesis of a
monotonic trend if any lower limit is greater than 0 and any upper limit
is less than 0. The sample median and standard error for each group can
be computed using the
[ci.median](https://dgbonett.github.io/statpsych/reference/ci.median.md)
function.

For more details, see Section 1.21 of Bonett (2021, Volume 2)

## Usage

``` r
test.mono.median.bs(alpha, m, se)
```

## Arguments

- alpha:

  alpha level for simultaneous 1-alpha confidence

- m:

  vector of estimated group medians

- se:

  vector of estimated group standard errors

## Value

Returns a matrix with the number of rows equal to the number of adjacent
pairwise comparisons. The columns are:

- Estimate - estimated median difference

- SE - standard error

- LL - one-sided lower limit of the confidence interval

- UL - one-sided upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
m <- c(12.86, 24.57, 36.29, 53.21)
se <- c(2.85, 2.99, 3.73, 3.88)
test.mono.median.bs(.05, m, se)
#>      Estimate       SE        LL         UL
#>  1 2   -11.71 4.130690 -21.59879 -1.8212115
#>  2 3   -11.72 4.780481 -23.16438 -0.2756247
#>  3 4   -16.92 5.382128 -29.80471 -4.0352947

# Should return:
#      Estimate       SE        LL         UL
#  1 2   -11.71 4.130690 -21.59879 -1.8212115
#  2 3   -11.72 4.780481 -23.16438 -0.2756247
#  3 4   -16.92 5.382128 -29.80471 -4.0352947

```
