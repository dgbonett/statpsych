# Computes confidence intervals of standardized effects in a 2x2 mixed design

Computes confidence intervals for the standardized AB interaction
effect, main effect of A, main effect of B, simple main effects of A,
and simple main effects of B in a 2x2 mixed factorial design where
Factor A is a within-subjects factor, and Factor B is a between-subjects
factor. Equality of population variances is not assumed. A square root
unweighted average variance standardizer is used.

## Usage

``` r
ci.2x2.stdmean.mixed(alpha, y11, y12, y21, y22)
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
y11 <- c(18, 19, 20, 17, 20, 16)
y12 <- c(19, 16, 16, 14, 16, 18)
y21 <- c(19, 18, 19, 20, 17, 16)
y22 <- c(16, 10, 12,  9, 13, 15)
ci.2x2.stdmean.mixed(.05, y11, y12, y21, y22)
#>          Estimate adj Estimate      SE      LL      UL
#> AB:       -1.9515      -1.8014 0.54074 -3.0114 -0.8917
#> A:         1.0606       1.0113 0.27977  0.5123  1.6090
#> B:         1.9091       1.7623 0.57589  0.7804  3.0378
#> A at b1:   0.0848       0.0759 0.46504 -0.8266  0.9963
#> A at b2:   2.0364       1.8214 0.29956  1.4493  2.6235
#> B at a1:   0.9333       0.8615 0.50364 -0.0538  1.9205
#> B at a2:   2.8849       2.6630 0.74772  1.4194  4.3504

# Should return:
#          Estimate adj Estimate      SE      LL      UL
# AB:       -1.9515      -1.8014 0.54074 -3.0114 -0.8917
# A:         1.0606       1.0113 0.27977  0.5123  1.6090
# B:         1.9091       1.7623 0.57589  0.7804  3.0378
# A at b1:   0.0848       0.0759 0.46504 -0.8266  0.9963
# A at b2:   2.0364       1.8214 0.29956  1.4493  2.6235
# B at a1:   0.9333       0.8615 0.50364 -0.0538  1.9205
# B at a2:   2.8849       2.6630 0.74772  1.4194  4.3504

```
