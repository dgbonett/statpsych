# Approximates the power of a correlation test for a planned sample size

Computes the approximate power of a test for a population Pearson or
partial correlation test for a planned sample size. Set s = 0 for a
Pearson correlation.

## Usage

``` r
power.cor(alpha, n, cor, h, s)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- n:

  planned sample size

- cor:

  planning value of correlation

- h:

  null hypothesis value of correlation

- s:

  number of control variables

## Value

Returns the approximate power of the test

## Examples

``` r
power.cor(.05, 80, .3, 0, 0)
#>      Power
#>  0.7751947

# Should return:
#     Power
# 0.7751947

```
