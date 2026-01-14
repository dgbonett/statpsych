# Computes confidence intervals of effects in a 2x2 within-subjects design for medians

Computes distribution-free confidence intervals for the AB interaction
effect, main effect of A, main effect of B, simple main effects of A,
and simple main effects of B in a 2x2 within-subjects factorial design.
The effects are defined in terms of medians rather than means. Tied
scores within each level combination are assumed to be rare.

## Usage

``` r
ci.2x2.median.ws(alpha, y11, y12, y21, y22)
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

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Price RM (2020). “Interval estimation for linear functions of
medians in within-subjects and mixed designs.” *British Journal of
Mathematical and Statistical Psychology*, **73**(2), 333–346. ISSN
0007-1102, [doi:10.1111/bmsp.12171](https://doi.org/10.1111/bmsp.12171)
.

## Examples

``` r
y11 <- c(222, 402, 333, 301, 284, 182, 281, 230, 290, 182, 133, 278)
y12 <- c(221, 371, 340, 288, 293, 150, 317, 211, 286, 161, 126, 234)
y21 <- c(219, 371, 314, 279, 284, 155, 278, 185, 296, 169, 118, 229)
y22 <- c(170, 332, 280, 273, 272, 160, 260, 204, 252, 153, 137, 223)
ci.2x2.median.ws(.05, y11, y12, y21, y22)
#>          Estimate        SE           LL       UL
#> AB:          3.50 21.050122 -37.75748155 44.75748
#> A:          24.25  9.603490   5.42750463 43.07250
#> B:          17.75  9.101881  -0.08935904 35.58936
#> A at b1:    26.00 11.813742   2.84549058 49.15451
#> A at b2:    22.50 16.323093  -9.49267494 54.49267
#> B at a1:    19.50 15.710347 -11.29171468 50.29171
#> B at a2:    16.00 11.850202  -7.22596953 39.22597

# Should return:
#          Estimate        SE           LL       UL
# AB:          3.50 21.050122 -37.75748155 44.75748
# A:          24.25  9.603490   5.42750463 43.07250
# B:          17.75  9.101881  -0.08935904 35.58936
# A at b1:    26.00 11.813742   2.84549058 49.15451
# A at b2:    22.50 16.323093  -9.49267494 54.49267
# B at a1:    19.50 15.710347 -11.29171468 50.29171
# B at a2:    16.00 11.850202  -7.22596953 39.22597

```
