# Computes tests and confidence intervals of effects in a 2x2 between-subjects design for medians

Computes distribution-free confidence intervals for the AB interaction
effect, main effect of A, main effect of B, simple main effects of A,
and simple main effects of B in a 2x2 between-subjects factorial design
with a quantitative response variable. The effects are defined in terms
of medians rather than means. Tied scores within each group are assumed
to be rare.

For more details, see Section 3.21 of Bonett (2021, Volume 1)

## Usage

``` r
ci.2x2.median.bs(alpha, y11, y12, y21, y22)
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

Bonett DG, Price RM (2002). “Statistical inference for a linear function
of medians: Confidence intervals, hypothesis testing, and sample size
requirements.” *Psychological Methods*, **7**(3), 370–383. ISSN
1939-1463,
[doi:10.1037/1082-989X.7.3.370](https://doi.org/10.1037/1082-989X.7.3.370)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
y11 <- c(19.2, 21.1, 14.4, 13.3, 19.8, 15.9, 18.0, 19.1, 16.2, 14.6)
y12 <- c(21.3, 27.0, 19.1, 21.5, 25.2, 24.1, 19.8, 19.7, 17.5, 16.0)
y21 <- c(16.5, 11.3, 10.3, 17.7, 13.8, 18.2, 12.8, 16.2, 6.1, 15.2)
y22 <- c(18.7, 17.3, 11.4, 12.4, 13.6, 13.8, 18.3, 15.0, 14.4, 11.9)
ci.2x2.median.bs(.05, y11, y12, y21, y22)
#>          Estimate       SE        LL         UL
#> AB:        -3.850 2.951019 -9.633891  1.9338914
#> A:          4.525 1.475510  1.633054  7.4169457
#> B:         -1.525 1.475510 -4.416946  1.3669457
#> A at b1:    2.600 1.992028 -1.304302  6.5043022
#> A at b2:    6.450 2.177232  2.182703 10.7172971
#> B at a1:   -3.450 2.045086 -7.458294  0.5582944
#> B at a2:    0.400 2.127472 -3.769769  4.5697694

# Should return:
#          Estimate       SE        LL         UL
# AB:        -3.850 2.951019 -9.633891  1.9338914
# A:          4.525 1.475510  1.633054  7.4169457
# B:         -1.525 1.475510 -4.416946  1.3669457
# A at b1:    2.600 1.992028 -1.304302  6.5043022
# A at b2:    6.450 2.177232  2.182703 10.7172971
# B at a1:   -3.450 2.045086 -7.458294  0.5582944
# B at a2:    0.400 2.127472 -3.769769  4.5697694

```
