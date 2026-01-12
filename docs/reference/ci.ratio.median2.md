# Confidence interval for a 2-group median ratio

Computes a distribution-free confidence interval for a ratio of
population medians of ratio-scale measurements in a 2-group design. Tied
scores within each group are assumed to be rare.

For more details, see Section 2.12 of Bonett (2021, Volume 1)

## Usage

``` r
ci.ratio.median2(alpha, y1, y2)
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

- Median1 - estimated median for group 1

- Median2 - estimated median for group 2

- Median1/Median2 - estimated ratio of medians

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Price RM (2020). “Confidence intervals for ratios of means
and medians.” *Journal of Educational and Behavioral Statistics*,
**45**(6), 750–770.
[doi:10.3102/1076998620934125](https://doi.org/10.3102/1076998620934125)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
y2 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29, 49, 42, 40)
y1 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
ci.ratio.median2(.05, y1, y2)
#>  Median1 Median2 Median1/Median2       LL       UL
#>       43      37        1.162162 0.927667 1.455933

# Should return:
# Median1 Median2 Median1/Median2       LL       UL
#      43      37        1.162162 0.927667 1.455933

```
