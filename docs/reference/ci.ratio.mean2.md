# Confidence interval for a 2-group mean ratio

Computes a confidence interval for a ratio of population means of
ratio-scale measurements in a 2-group design. Equality of variances is
not assumed.

## Usage

``` r
ci.ratio.mean2(alpha, y1, y2)
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

- Mean1 - estimated mean for group 1

- Mean2 - estimated mean for group 2

- Mean1/Mean2- estimated mean ratio

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## Details

For more details, see Section 2.5 of Bonett (2021, Volume 1)

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
ci.ratio.mean2(.05, y1, y2)
#>  Mean1    Mean2 Mean1/Mean2        LL      UL
#>   41.5 36.38462    1.140592 0.9837277 1.32247

# Should return:
#
# Mean1    Mean2 Mean1/Mean2        LL       UL
#  41.5 36.38462    1.140592 0.9897482 1.314425

```
