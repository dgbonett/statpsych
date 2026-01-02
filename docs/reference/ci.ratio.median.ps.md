# Confidence interval for a paired-samples median ratio

Computes a distribution-free confidence interval for a ratio of
population medians in a paired-samples design. Ratio-scale measurements
are assumed. Tied scores within each measurement are assumed to be rare.

For more details, see Section 4.23 of Bonett (2021, Volume 1)

## Usage

``` r
ci.ratio.median.ps(alpha, y1, y2)
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

- Median1/Median2 - estimated ratio of medians

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
ci.ratio.median.ps(.05, y1, y2)
#>  Median1 Median2 Median1/Median2       LL       UL
#>   81.835   67.73        1.208253 1.069251 1.365326

# Should return:
# Median1  Median2  Median1/Median2        LL        UL
#  81.835    67.73         1.208253  1.069251  1.365326

```
