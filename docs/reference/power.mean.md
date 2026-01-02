# Approximates the power of a one-sample t-test for a planned sample size

Computes the approximate power of a one-sample t-test for a planned
sample size. For a conservatively low power approximation, set the
variance planning value to the largest value within its plausible range,
and set the effect size to a minimally interesting value.

## Usage

``` r
power.mean(alpha, n, var, es)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- n:

  planned sample size

- var:

  planning value of response variable variance

- es:

  planning value of mean minus null hypothesis value

## Value

Returns the approximate power of the test

## Examples

``` r
power.mean(.05, 15, 80.5, 7)
#>      Power
#>  0.8021669

# Should return:
#     Power
# 0.8021669

```
