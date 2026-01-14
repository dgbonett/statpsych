# Confidence interval for two kappa reliability coefficients

Computes confidence intervals for the intraclass kappa coefficient and
Cohen's kappa coefficient with two dichotomous ratings. The G-index of
agreement (see ci.agree) is arguably a better measure of agreement.

For more details, see Section 3.5 of Bonett (2021, Volume 3)

## Usage

``` r
ci.kappa(alpha, f00, f01, f10, f11)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f00:

  number of objects rated 0 by both Rater 1 and Rater 2

- f01:

  number of objects rated 0 by Rater 1 and 1 by Rater 2

- f10:

  number of objects rated 1 by Rater 1 and 0 by Rater 2

- f11:

  number of objects rated 1 by both Rater 1 and Rater 2

## Value

Returns a 2-row matrix. The results in row 1 are for the intraclass
kappa. The results in row 2 are for Cohen's kappa. The columns are:

- Estimate - estimate of interrater reliability

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Fleiss JL, Paik MC (2003). *Statistical Methods for Rates and
Proportions*, 3rd edition. Wiley.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.kappa(.05, 31, 12, 4, 58)
#>              Estimate     SE    LL    UL
#> IC kappa:       0.674 0.0748 0.527 0.821
#> Cohen kappa:    0.676 0.0734 0.532 0.820

# Should return:
#              Estimate     SE    LL    UL
# IC kappa:       0.674 0.0748 0.527 0.821
# Cohen kappa:    0.676 0.0734 0.532 0.820

```
