# Computes confidence intervals of standardized effects in a 2x2 between-subjects design

Computes confidence intervals for standardized AB interaction effect,
main effect of A, main effect of B, simple main effects of A, and simple
main effects of B in a 2x2 between-subjects factorial design with a
quantitative response variable. Equality of population variances is not
assumed. A square root unweighted average variance standardizer is used,
which is the recommended standardizer when both factors are treatment
factors.

## Usage

``` r
ci.2x2.stdmean.bs(alpha, y11, y12, y21, y22)
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

- Estimate - estimate of standardized effect

- adj Estimate - bias adjusted estimate of standardized effect

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
y11 <- c(14, 15, 11, 7, 16, 12, 15, 16, 10, 9)
y12 <- c(18, 24, 14, 18, 22, 21, 16, 17, 14, 13)
y21 <- c(16, 11, 10, 17, 13, 18, 12, 16, 6, 15)
y22 <- c(18, 17, 11, 9, 9, 13, 18, 15, 14, 11)
ci.2x2.stdmean.bs(.05, y11, y12, y21, y22)
#>          Estimate adj Estimate      SE      LL      UL
#> AB:       -1.4498      -1.4194 0.68852 -2.7992 -0.1003
#> A:         0.4690       0.4592 0.33795 -0.1933  1.1314
#> B:        -0.7533      -0.7375 0.34512 -1.4297 -0.0769
#> A at b1:  -0.2558      -0.2505 0.46402 -1.1653  0.6536
#> A at b2:   1.1939       1.1689 0.50014  0.2137  2.1742
#> B at a1:  -1.4782      -1.4472 0.49284 -2.4441 -0.5122
#> B at a2:  -0.0284      -0.0278 0.48204 -0.9732  0.9163

# Should return:
#          Estimate adj Estimate      SE      LL      UL
# AB:       -1.4498      -1.4194 0.68852 -2.7992 -0.1003
# A:         0.4690       0.4592 0.33795 -0.1933  1.1314
# B:        -0.7533      -0.7375 0.34512 -1.4297 -0.0769
# A at b1:  -0.2558      -0.2505 0.46402 -1.1653  0.6536
# A at b2:   1.1939       1.1689 0.50014  0.2137  2.1742
# B at a1:  -1.4782      -1.4472 0.49284 -2.4441 -0.5122
# B at a2:  -0.0284      -0.0278 0.48204 -0.9732  0.9163

```
