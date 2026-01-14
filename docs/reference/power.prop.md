# Approximates the power of a 1-group proportion test for a planned sample size

Computes the approximate power of a one-sample proportion test for a
planned sample size. Set the proportion planning value to .5 for a
conservatively low power estimate. The value of the effect size need not
be based on the proportion planning value.

## Usage

``` r
power.prop(alpha, n, p, es)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- n:

  planned sample size

- p:

  planning value of proportion

- es:

  planning value of proportion minus null hypothesis value

## Value

Returns the approximate power of the test

## Examples

``` r
power.prop(.05, 40, .5, .2)
#>      Power
#>  0.7156166

# Should return:
#     Power
# 0.7156044

```
