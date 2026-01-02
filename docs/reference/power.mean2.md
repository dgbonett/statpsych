# Approximates the power of a two-sample t-test for planned sample sizes

Computes the approximate power of a two-sample t-test for planned sample
sizes. For a conservatively low power approximation, set the variance
planning values to the largest values within their plausible ranges, and
set the effect size to a minimally interesting value. The within-group
variances can be unequal across groups and a Satterthwaite degree of
freedom adjustment is used to improve the accuracy of the power
approximation.

## Usage

``` r
power.mean2(alpha, n1, n2, var1, var2, es)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- n1:

  planned sample size for group 1

- n2:

  planned sample size for group 2

- var1:

  planning value of within-group variance for group 1

- var2:

  planning value of within-group variance for group 2

- es:

  planning value of mean difference

## Value

Returns the approximate power of the test

## Examples

``` r
power.mean2(.05, 25, 25, 5.0, 6.0, 2)
#>      Power
#>  0.8398417

# Should return:
#     Power
# 0.8398417

```
