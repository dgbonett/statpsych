# Approximates the power of a test for equal correlations in a 2-group design for planned sample sizes

Computes the approximate power of a test for equal population Pearson or
partial correlations in a 2-group design for planned sample sizes. Set s
= 0 for Pearson correlations.

## Usage

``` r
power.cor2(alpha, n1, n2, cor1, cor2, s)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- n1:

  planned sample size for group 1

- n2:

  planned sample size for group 2

- cor1:

  correlation planning value for group 1

- cor2:

  correlation planning value for group 2

- s:

  number of control variables

## Value

Returns the approximate power of the test

## Examples

``` r
power.cor2(.05, 200, 200, .4, .2, 0)
#>      Power
#>  0.5919682

# Should return:
#     Power
# 0.5919682

```
