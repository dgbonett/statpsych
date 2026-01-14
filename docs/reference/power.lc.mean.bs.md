# Approximates the power of a test for a linear contrast of means for planned sample sizes in a between-subjects design

Computes the approximate power of a test for a linear contrast of
population means for planned sample sizes in a between-subject design.
The groups can be the factor levels of a single factor design or the
combinations of factors in a factorial design. For a conservatively low
power approximation, set the variance planning values to the largest
values within their plausible ranges, and set the effect size to a
minimally interesting value. The within-group variances can be unequal
across groups and a Satterthwaite degree of freedom adjustment is used
to improve the accuracy of the power approximation.

## Usage

``` r
power.lc.mean.bs(alpha, n, var, es, v)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- n:

  vector of planned sample sizes

- var:

  vector of within-group variance planning values

- es:

  planning value of linear contrast of means

- v:

  vector of contrast coefficients

## Value

Returns the approximate power of the test

## Examples

``` r
n <- c(20, 20, 20, 20)
var <- c(70, 70, 80, 80)
v <- c(.5, .5, -.5, -.5)
power.lc.mean.bs(.05, n, var, 5, v)
#>      Power
#>  0.7221171

# Should return:
#     Power
# 0.7221171

```
