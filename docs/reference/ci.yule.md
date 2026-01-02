# Confidence intervals for generalized Yule coefficients

Computes confidence intervals for four generalized Yule measures of
association (Yule Q, Yule Y, Digby H, and Bonett-Price Y\*) using a
transformation of a confidence interval for an odds ratio with .5 added
to each cell frequency. This function requires the frequency counts from
a 2 x 2 contingency table for two dichotomous variables. Digby H is
sometimes used as a crude approximation to the tetrachoric correlation.
Yule Y is equal to the phi coefficient only when all marginal
frequencies are equal. Bonett-Price Y\* is a better approximation to the
phi coefficient when the marginal frequencies are not equal.

For more details, see Section 3.4 of Bonett (2021, Volume 3)

## Usage

``` r
ci.yule(alpha, f00, f01, f10, f11)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f00:

  number of participants with y = 0 and x = 0

- f01:

  number of participants with y = 0 and x = 1

- f10:

  number of participants with y = 1 and x = 0

- f11:

  number of participants with y = 1 and x = 1

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of generalized Yule coefficient

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Price RM (2007). “Statistical inference for generalized Yule
coefficients in 2x2 contingency tables.” *Sociological Methods &
Research*, **35**(3), 429–446. ISSN 0049-1241,
[doi:10.1177/0049124106292358](https://doi.org/10.1177/0049124106292358)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.yule(.05, 229, 28, 96, 24)
#>     Estimate     SE    LL    UL
#> Q:     0.343 0.1328 0.062 0.573
#> Y:     0.177 0.0729 0.031 0.315
#> H:     0.262 0.1051 0.047 0.454
#> Y*:    0.131 0.0546 0.023 0.236

# Should return:
#    Estimate     SE    LL    UL
# Q:    0.343 0.1328 0.062 0.573
# Y:    0.177 0.0729 0.031 0.315
# H:    0.262 0.1051 0.047 0.454
# Y*:   0.131 0.0546 0.023 0.236

```
