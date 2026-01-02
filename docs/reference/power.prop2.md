# Approximates the power of a 2-group proportion test for planned sample sizes

Computes the approximate power for a test of equal population
proportions in a 2-group design for the planned sample sizes. This
function requires planning values for both proportions. Set the
proportion planning values to .5 for a conservatively low power
estimate. The planning value for the proportion difference could be set
to the difference of the two proportion planning values or it could be
set to a minimally interesting effect size.

## Usage

``` r
power.prop2(alpha, n1, n2, p1, p2, es)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- n1:

  planned sample size for group 1

- n2:

  planned sample size for group 2

- p1:

  planning value of proportion for group 1

- p2:

  planning value of proportion for group 2

- es:

  planning value of proportion difference

## Value

Returns the approximate power of the test

## Examples

``` r
power.prop2(.05, 60, 40, .5, .5, .2)
#>      Power
#>  0.4998959

# Should return:
#     Power
# 0.4998959

```
