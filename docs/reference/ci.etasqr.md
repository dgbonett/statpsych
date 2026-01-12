# Confidence interval for eta-squared

Computes a confidence interval for a population eta-squared, partial
eta-squared, or generalized eta-squared in a fixed-factor
between-subjects design. An approximate bias adjusted estimate is
computed, and an approximate standard error is recovered from the
confidence interval.

For more details, see Section 3.7 of Bonett (2021, Volume 1)

## Usage

``` r
ci.etasqr(alpha, etasqr, df1, df2)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- etasqr:

  estimated eta-squared

- df1:

  degrees of freedom for effect

- df2:

  error degrees of freedom

## Value

Returns a 1-row matrix. The columns are:

- Eta-squared - eta-squared (from input)

- adj Eta-squared - bias adjusted eta-squared estimate

- SE - recovered standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.etasqr(.05, .15, 2, 57)
#>  Eta-squared adj Eta-squared      SE     LL     UL
#>         0.15          0.1202 0.07414 0.0102 0.3008

# Should return:
# Eta-squared  adj Eta-squared       SE      LL      UL
#        0.15           0.1202  0.07414  0.0102  0.3008
 
```
