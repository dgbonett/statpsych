# Computes confidence intervals in a 2x2 mixed design for medians

Computes distribution-free confidence intervals for the AB interaction
effect, main effect of A, main effect of B, simple main effects of A,
and simple main effects of B in a 2x2 mixed factorial design where
Factor A is the within-subjects factor and Factor B is the between
subjects factor. Tied scores within each group and within each
within-subjects level are assumed to be rare.

## Usage

``` r
ci.2x2.median.mixed(alpha, y11, y12, y21, y22)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y11:

  vector of scores at level 1 of A and level 1 of B (group 1)

- y12:

  vector of scores at level 1 of A and level 2 of B (group 2)

- y21:

  vector of scores at level 2 of A and level 1 of B (group 1)

- y22:

  vector of scores at level 2 of A and level 2 of B (group 2)

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
y11 <- c(18.3, 19.5, 20.1, 17.4, 20.5, 16.1)
y12 <- c(19.2, 16.4, 16.5, 14.0, 16.9, 18.3)
y21 <- c(19.1, 18.4, 19.8, 20.0, 17.2, 16.8)
y22 <- c(16.5, 10.2, 12.7,  9.9, 13.5, 15.0)
ci.2x2.median.mixed(.05, y11, y12, y21, y22)
#>          Estimate        SE         LL         UL
#> AB:        -3.450 1.6317863 -6.6482423 -0.2517577
#> A:          3.925 0.8158931  2.3258788  5.5241212
#> B:          1.875 1.4262367 -0.9203726  4.6703726
#> A at b1:    0.150 1.4243192 -2.6416144  2.9416144
#> A at b2:    3.600 0.7962670  2.0393454  5.1606546
#> B at a1:    2.200 1.5812792 -0.8992503  5.2992503
#> B at a2:    5.650 1.7027101  2.3127496  8.9872504

# Should return:
#          Estimate        SE         LL         UL
# AB:        -3.450 1.6317863 -6.6482423 -0.2517577
# A:          1.875 0.8158931  0.2758788  3.4741212
# B:          3.925 1.4262367  1.1296274  6.7203726
# A at b1:    0.150 1.4243192 -2.6416144  2.9416144
# A at b2:    3.600 0.7962670  2.0393454  5.1606546
# B at a1:    2.200 1.5812792 -0.8992503  5.2992503
# B at a2:    5.650 1.7027101  2.3127496  8.9872504

```
