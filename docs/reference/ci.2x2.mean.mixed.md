# Computes tests and confidence intervals of effects in a 2x2 mixed design for means

Computes confidence intervals and tests for the AB interaction effect,
main effect of A, main effect of B, simple main effects of A, and simple
main effects of B in a 2x2 mixed factorial design with a quantitative
response variable where Factor A is a within-subjects factor and Factor
B is a between-subjects factor. A Satterthwaite adjustment to the
degrees of freedom is used and equality of population variances is not
assumed.

For more details, see Section 4.16 of Bonett (2021, Volume 1)

## Usage

``` r
ci.2x2.mean.mixed(alpha, y11, y12, y21, y22)
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

- t - t test statistic

- df - degrees of freedom

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
y11 <- c(18, 19, 20, 17, 20, 16)
y12 <- c(19, 16, 16, 14, 16, 18)
y21 <- c(19, 18, 19, 20, 17, 16)
y22 <- c(16, 10, 12,  9, 13, 15)
ci.2x2.mean.mixed(.05, y11, y12, y21, y22)
#>            Estimate        SE       t   df       p         LL        UL
#> AB:      -3.8333333 0.9803627 -3.9101 8.35 0.00412 -6.0776652 -1.589001
#> A:        2.0833333 0.4901814  4.2501 8.35 0.00254  0.9611674  3.205499
#> B:        3.7500000 1.0226599  3.6669 7.60 0.00693  1.3699617  6.130038
#> A at b1:  0.1666667 0.8333333  0.2000 5.00 0.84936 -1.9754849  2.308818
#> A at b2:  4.0000000 0.5163978  7.7460 5.00 0.00057  2.6725572  5.327443
#> B at a1:  1.8333333 0.9803627  1.8701 9.94 0.09117 -0.3528396  4.019506
#> B at a2:  5.6666667 1.2692955  4.4644 7.67 0.00233  2.7176000  8.615733

# Should return:
#            Estimate        SE       t   df       p         LL        UL
# AB:      -3.8333333 0.9803627 -3.9101 8.35 0.00412 -6.0778198 -1.588847
# A:        2.0833333 0.4901814  4.2501 8.35 0.00254  0.9610901  3.205577
# B:        3.7500000 1.0226599  3.6669 7.60 0.00693  1.3700362  6.129964
# A at b1:  0.1666667 0.8333333  0.2000 5.00 0.84936 -1.9754849  2.308818
# A at b2:  4.0000000 0.5163978  7.7460 5.00 0.00057  2.6725572  5.327443
# B at a1:  1.8333333 0.9803627  1.8701 9.94 0.09117 -0.3527241  4.019391
# B at a2:  5.6666667 1.2692955  4.4644 7.67 0.00233  2.7173445  8.615989

```
