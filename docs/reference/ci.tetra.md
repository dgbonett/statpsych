# Confidence interval for a tetrachoric correlation

Computes a confidence interval for an approximation to the tetrachoric
correlation. This function requires the frequency counts from a 2 x 2
contingency table for two dichotomous variables. This measure of
association assumes both of the dichotomous variables are artificially
dichotomous. An approximate standard error is recovered from the
confidence interval.

For more details, see Section 3.4 of Bonett (2021, Volume 3)

## Usage

``` r
ci.tetra(alpha, f00, f01, f10, f11)
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

- Estimate - estimate of tetrachoric approximation

- SE - recovered standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

Bonett DG, Price RM (2005). “Inferential methods for the tetrachoric
correlation coefficient.” *Journal of Educational and Behavioral
Statistics*, **30**(2), 213–225. ISSN 1076-9986,
[doi:10.3102/10769986030002213](https://doi.org/10.3102/10769986030002213)
.

## Examples

``` r
ci.tetra(.05, 86, 16, 7, 93)
#>  Estimate     SE    LL    UL
#>     0.938 0.0268 0.868 0.973

# Should return:
# Estimate     SE    LL    UL
#    0.938 0.0268 0.868 0.973

```
