# Confidence interval for a paired-samples proportion difference

Computes an adjusted Wald confidence interval for a difference of
population proportions in a paired-samples design. This function
requires the frequency counts from a 2 x 2 contingency table for two
repeated dichotomous measurements.

For more details, see Section 3.2 of Bonett (2021, Volume 3)

## Usage

``` r
ci.prop.ps(alpha, f00, f01, f10, f11)
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

- Estimate - adjusted estimate of proportion difference

- SE - adjusted standard error

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Bonett DG, Price RM (2012). “Adjusted wald confidence interval for a
difference of binomial proportions based on paired data.” *Journal of
Educational and Behavioral Statistics*, **37**(4), 479–488. ISSN
1076-9986,
[doi:10.3102/1076998611411915](https://doi.org/10.3102/1076998611411915)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.prop.ps(.05, 12, 4, 26, 6)
#>  Estimate         SE        LL        UL
#>      0.44 0.09448809 0.2548067 0.6251933

# Should return:
# Estimate         SE        LL        UL
#     0.44 0.09448809 0.2548067 0.6251933

```
