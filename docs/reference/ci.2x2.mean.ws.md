# Computes tests and confidence intervals of effects in a 2x2 within-subjects design for means

Computes confidence intervals and tests for the AB interaction effect,
main effect of A, main effect of B, simple main effects of A, and simple
main effects of B in a 2x2 within-subjects factorial design with a
quantitative response variable.

For more details, see Section 4.14 of Bonett (2021, Volume 1)

## Usage

``` r
ci.2x2.mean.ws(alpha, y11, y12, y21, y22)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y11:

  vector of scores at level 1 of A and level 1 of B

- y12:

  vector of scores at level 1 of A and level 2 of B

- y21:

  vector of scores at level 2 of A and level 1 of B

- y22:

  vector of scores at level 2 of A and level 2 of B

## Value

Returns a 7-row matrix (one row per effect). The columns are:

- Estimate - estimate of effect

- SE - standard error

- t - t test statistic

- df - degrees of freedom

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
y11 <- c(1,2,3,4,5,7,7)
y12 <- c(1,0,2,4,3,8,7)
y21 <- c(4,5,6,7,8,9,8)
y22 <- c(5,6,8,7,8,9,9)
ci.2x2.mean.ws(.05, y11, y12, y21, y22)
#>             Estimate        SE       t df       p          LL          UL
#> AB:       1.28571429 0.5654449  2.2738  6 0.06334 -0.09787945  2.66930802
#> A:       -3.21428571 0.4862042 -6.6110  6 0.00058 -4.40398462 -2.02458681
#> B:       -0.07142857 0.2296107 -0.3111  6 0.76625 -0.63326579  0.49040865
#> A at b1: -2.57142857 0.2973809 -8.6469  6 0.00013 -3.29909331 -1.84376383
#> A at b2: -3.85714286 0.7377111 -5.2285  6 0.00196 -5.66225692 -2.05202879
#> B at a1:  0.57142857 0.4285714  1.3333  6 0.23082 -0.47724794  1.62010508
#> B at a2: -0.71428571 0.2857143 -2.5000  6 0.04653 -1.41340339 -0.01516804

# Should return:
#             Estimate        SE       t df       p          LL          UL
# AB:       1.28571429 0.5654449  2.2738  6 0.06334 -0.09787945  2.66930802
# A:       -3.21428571 0.4862042 -6.6110  6 0.00058 -4.40398462 -2.02458681
# B:       -0.07142857 0.2296107 -0.3111  6 0.76626 -0.63326579  0.49040865
# A at b1: -2.57142857 0.2973809 -8.6469  6 0.00013 -3.29909331 -1.84376383
# A at b2: -3.85714286 0.7377111 -5.2285  6 0.00196 -5.66225692 -2.05202879
# B at a1:  0.57142857 0.4285714  1.3333  6 0.23081 -0.47724794  1.62010508
# B at a2: -0.71428571 0.2857143 -2.5000  6 0.04653 -1.41340339 -0.01516804

```
