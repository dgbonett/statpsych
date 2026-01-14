# Approximates the power of a paired-samples test of equal proportions for a planned sample size

Computes the approximate power of a test for equal population
proportions in a paired-samples design (the McNemar test). This function
requires planning values for both proportions and a phi coefficient that
describes the correlation between the two dichotomous measurements. The
proportion planning values can set to .5 for a conservatively low power
estimate. The planning value for the proportion difference (effect size)
could be set to the difference of the two proportion planning values or
it could be set to a minimally interesting effect size. Set the phi
correlation planning value to the smallest value within a plausible
range for a conservatively low power estimate.

## Usage

``` r
power.prop.ps(alpha, n, p1, p2, phi, es)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- n:

  planned sample size

- p1:

  planning value of proportion for measurement 1

- p2:

  planning value of proportion for measurement 2

- phi:

  planning value of phi correlation

- es:

  planning value of proportion difference

## Value

Returns the approximate power of the test

## Examples

``` r
power.prop.ps(.05, 45, .5, .5, .4, .2)
#>      Power
#>  0.6877704

# Should return:
#     Power
# 0.6877704

```
