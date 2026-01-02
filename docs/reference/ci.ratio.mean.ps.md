# Confidence interval for a paired-samples mean ratio

Compute a confidence interval for a ratio of population means of
ratio-scale measurements in a paired-samples design. Equality of
variances is not assumed.

For more details, see Section 4.4 of Bonett (2021, Volume 1)

## Usage

``` r
ci.ratio.mean.ps(alpha, y1, y2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y1:

  vector of measurement 1 scores

- y2:

  vector of measurement 2 scores (paired with y1)

## Value

Returns a 1-row matrix. The columns are:

- Mean1 - estimated mean for measurement 1

- Mean2 - estimated mean for measurement 2

- Mean1/Mean2 - estimate of mean ratio

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Price RM (2020). “Confidence intervals for ratios of means
and medians.” *Journal of Educational and Behavioral Statistics*,
**45**(6), 750–770.
[doi:10.3102/1076998620934125](https://doi.org/10.3102/1076998620934125)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
y1 = c(76.41, 66.91, 81.06, 74.78, 83.76, 89.31, 78.78, 87.06, 82.61, 76.74, 88.33, 86.18)
y2 = c(59.85, 60.64, 84.86, 68.16, 71.53, 86.18, 67.30, 65.46, 83.50, 66.76, 88.37, 65.02)
ci.ratio.mean.ps(.05, y1, y2)
#>     Mean1   Mean2 Mean1/Mean2       LL       UL
#>  80.99417 72.3025    1.120213 1.040747 1.205745

# Should return:
#    Mean1   Mean2 Mean1/Mean2       LL       UL
# 80.99417 72.3025    1.120213 1.040747 1.205745

```
