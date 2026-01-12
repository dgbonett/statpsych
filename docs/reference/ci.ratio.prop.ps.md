# Confidence interval for a paired-samples proportion ratio

Computes a confidence interval for a ratio of population proportions in
a paired-samples design. This function requires the frequency counts
from a 2 x 2 contingency table for two repeated dichotomous
measurements.

For more details, see Section 3.2 of Bonett (2021, Volume 3)

## Usage

``` r
ci.ratio.prop.ps(alpha, f00, f01, f10, f11)
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

- Estimate - estimate of proportion ratio

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Price RM (2006). “Confidence intervals for a ratio of
binomial proportions based on paired data.” *Statistics in Medicine*,
**25**(17), 3039–3047. ISSN 0277-6715,
[doi:10.1002/sim.2440](https://doi.org/10.1002/sim.2440) .

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.ratio.prop.ps(.05, 12, 4, 26, 6)
#>  Estimate       LL       UL
#>       3.2 1.766544 5.796628

# Should return:
# Estimate       LL       UL
#      3.2 1.766544 5.796628

```
