# Computes tests and confidence intervals of effects in a 2x2 between-subjects design for means

Computes confidence intervals and tests for the AB interaction effect,
main effect of A, main effect of B, simple main effects of A, and simple
main effects of B in a 2x2 between-subjects factorial design with a
quantitative response variable. A Satterthwaite adjustment to the
degrees of freedom is used and equality of population variances is not
assumed.

For more details, see Sections 3.8 and 3.9 of Bonett (2021, Volume 1)

## Usage

``` r
ci.2x2.mean.bs(alpha, y11, y12, y21, y22)
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

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
y11 <- c(14, 15, 11, 7, 16, 12, 15, 16, 10, 9)
y12 <- c(18, 24, 14, 18, 22, 21, 16, 17, 14, 13)
y21 <- c(16, 11, 10, 17, 13, 18, 12, 16, 6, 15)
y22 <- c(18, 17, 11, 9, 9, 13, 18, 15, 14, 11)
ci.2x2.mean.bs(.05, y11, y12, y21, y22)
#>          Estimate       SE       t    df       p         LL         UL
#> AB:         -5.10 2.224860 -2.2923 35.48 0.02793 -9.6145217 -0.5854783
#> A:           1.65 1.112430  1.4832 35.48 0.14685 -0.6072608  3.9072608
#> B:          -2.65 1.112430 -2.3822 35.48 0.02270 -4.9072608 -0.3927392
#> A at b1:    -0.90 1.545244 -0.5824 17.56 0.56770 -4.1522770  2.3522770
#> A at b2:     4.20 1.600694  2.6239 17.94 0.01724  0.8362597  7.5637403
#> B at a1:    -5.20 1.536952 -3.3833 17.61 0.00339 -8.4341504 -1.9658496
#> B at a2:    -0.10 1.608657 -0.0622 17.92 0.95109 -3.4807450  3.2807450

# Should return:
#          Estimate       SE       t    df       p         LL         UL
# AB:         -5.10 2.224860 -2.2923 35.48 0.02793 -9.6145264 -0.5854736
# A:           1.65 1.112430  1.4832 35.48 0.14684 -0.6072632  3.9072632
# B:          -2.65 1.112430 -2.3822 35.48 0.02270 -4.9072632 -0.3927368
# A at b1:    -0.90 1.545244 -0.5824 17.56 0.56768 -4.1522367  2.3522367
# A at b2:     4.20 1.600694  2.6239 17.94 0.01725  0.8362274  7.5637726
# B at a1:    -5.20 1.536952 -3.3833 17.61 0.00339 -8.4341379 -1.9658621
# B at a2:    -0.10 1.608657 -0.0622 17.92 0.95112 -3.4807927  3.2807927

```
