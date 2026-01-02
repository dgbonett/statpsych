# Computes confidence intervals of standardized effects in a 2x2 within-subjects design

Computes confidence intervals for standardized AB interaction effect,
main effect of A, main effect of B, simple main effects of A, and simple
main effects of B in a 2x2 within-subjects factorial design. Equality of
population variances is not assumed. A square root unweighted average
variance standardizer is used.

## Usage

``` r
ci.2x2.stdmean.ws(alpha, y11, y12, y21, y22)
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

- Estimate - estimated standardized effect

- adj Estimate - bias adjusted standardized effect estimate

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2008). “Confidence intervals for standardized linear
contrasts of means.” *Psychological Methods*, **13**(2), 99–109. ISSN
1939-1463,
[doi:10.1037/1082-989X.13.2.99](https://doi.org/10.1037/1082-989X.13.2.99)
.

## Examples

``` r
y11 <- c(21, 39, 32, 29, 27, 17, 27, 21, 28, 17, 12, 27)
y12 <- c(20, 36, 33, 27, 28, 14, 30, 20, 27, 15, 11, 22)
y21 <- c(21, 36, 30, 27, 28, 15, 27, 18, 29, 16, 11, 22)
y22 <- c(18, 34, 29, 28, 28, 17, 27, 21, 26, 16, 14, 23)
ci.2x2.stdmean.ws(.05, y11, y12, y21, y22)
#>          Estimate adj Estimate      SE      LL     UL
#> AB:        0.1725       0.1645 0.13655 -0.0951 0.4401
#> A:         0.1092       0.1042 0.05753 -0.0035 0.2220
#> B:         0.0747       0.0713 0.05921 -0.0413 0.1908
#> A at b1:   0.1955       0.1864 0.08461  0.0297 0.3613
#> A at b2:   0.0230       0.0219 0.09372 -0.1607 0.2067
#> B at a1:   0.1610       0.1535 0.09457 -0.0244 0.3463
#> B at a2:  -0.0115      -0.0110 0.08596 -0.1800 0.1570

# Should return:
#          Estimate adj Estimate      SE      LL     UL
# AB:        0.1725       0.1645 0.13655 -0.0951 0.4401
# A:         0.1092       0.1042 0.05753 -0.0035 0.2220
# B:         0.0747       0.0713 0.05921 -0.0413 0.1908
# A at b1:   0.1955       0.1864 0.08461  0.0297 0.3613
# A at b2:   0.0230       0.0219 0.09372 -0.1607 0.2067
# B at a1:   0.1610       0.1535 0.09457 -0.0244 0.3463
# B at a2:  -0.0115      -0.0110 0.08596 -0.1800 0.1570

```
