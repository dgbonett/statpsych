# Approximates the power of a paired-samples t-test for a planned sample size

Computes the approximate power of a paired-samples t-test for a planned
sample size. For a conservatively low power approximation, set the
variance planning values to the largest values within their plausible
ranges, set the correlation planning value to the smallest value within
its plausible range, and set the effect size to a minimally interesting
value. The variances of the two measurements can be unequal.

## Usage

``` r
power.mean.ps(alpha, n, var1, var2, es, cor)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- n:

  planned sample size

- var1:

  planning value of measurement 1 variance

- var2:

  planning value of measurement 2 variance

- es:

  planning value of mean difference

- cor:

  planning value of correlation between measurements

## Value

Returns the approximate power of the test

## Examples

``` r
power.mean.ps(.05, 20, 10.0, 12.0, 2, .7)
#>      Power
#>  0.9074354

# Should return:
#     Power
# 0.9074354

```
