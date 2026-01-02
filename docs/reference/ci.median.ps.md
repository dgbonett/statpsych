# Confidence interval for a paired-samples median difference

Computes a distribution-free confidence interval for a difference of
population medians in a paired-samples design. This function also
computes the standard error of each median and the covariance between
the two estimated medians. Tied scores within each measurement are
assumed to be rare.

For more details, see Section 4.23 of Bonett (2021, Volume 1)

## Usage

``` r
ci.median.ps(alpha, y1, y2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y1:

  vector of scores for measurement 1

- y2:

  vector of scores for measurement 2 (paired with y1)

## Value

Returns a 1-row matrix. The columns are:

- Median1 - estimated median for measurement 1

- Median2 - estimated median for measurement 2

- Median1-Median2 - estimated difference of medians

- SE1 - standard error of median 1

- SE2 - standard error of median 2

- COV - covariance of the two estimated medians

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Price RM (2020). “Interval estimation for linear functions of
medians in within-subjects and mixed designs.” *British Journal of
Mathematical and Statistical Psychology*, **73**(2), 333–346. ISSN
0007-1102, [doi:10.1111/bmsp.12171](https://doi.org/10.1111/bmsp.12171)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
y1 <- c(21.1, 4.9, 9.2, 12.4, 35.8, 18.1, 10.7, 22.9, 24.0, 1.2, 6.1, 8.3, 13.1, 16.2)
y2 <- c(67.0, 28.1, 30.9, 28.6, 52.0, 40.8, 25.8, 37.4, 44.9, 10.3, 14.9, 20.2, 28.8, 40.6)
ci.median.ps(.05, y1, y2)
#>  Median1 Median2 Median1-Median2       SE        LL        UL      SE1      SE2
#>    12.75   29.85           -17.1 3.704248 -24.36019 -9.839807 3.379695 4.968956
#>      COV
#>  11.1957

# Should return:
#  Median1 Median2 Median1-Median2       SE        LL        UL
#    12.75   29.85           -17.1 3.704248 -24.36019 -9.839807
#       SE1      SE2     COV
#  3.379695 4.968956 11.1957

```
