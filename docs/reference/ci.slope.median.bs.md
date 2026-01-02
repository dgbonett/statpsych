# Confidence interval for the slope of medians in a one-factor experimental design with a quantitative between-subjects factor

Computes a distribution-free test and confidence interval for the slope
of medians in a one-factor experimental design with a quantitative
between-subjects factor using sample group medians and standard errors
as input. The sample median and standard error for each group can be
computed using the
[ci.median](https://dgbonett.github.io/statpsych/reference/ci.median.md)
function.

## Usage

``` r
ci.slope.median.bs(alpha, m, se, x)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m:

  vector of sample median

- se:

  vector of standard errors

- x:

  vector of quantitative factor values

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated slope

- SE - standard error

- z - z test statistic

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## Examples

``` r
m <- c(33.5, 37.9, 38.0, 44.1)
se <- c(0.84, 0.94, 1.65, 2.98)
x <- c(5, 10, 20, 30)
ci.slope.median.bs(.05, m, se, x)
#>   Estimate        SE      z       p        LL        UL
#>  0.3664407 0.1163593 3.1492 0.00164 0.1383806 0.5945008

# Should return:
#   Estimate        SE      z       p        LL        UL
#  0.3664407 0.1163593 3.1492 0.00164 0.1383806 0.5945008

```
