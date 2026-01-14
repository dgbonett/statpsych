# Confidence intervals for a paired-samples standardized mean difference

Computes confidence intervals for a population standardized mean
difference in a paired-samples design. A square root unweighted variance
standardizer and single measurement standard deviation standardizers are
used. Equality of variances is not assumed.

For more details, see Section 4.3 of Bonett (2021, Volume 1)

## Usage

``` r
ci.stdmean.ps(alpha, m1, m2, sd1, sd2, cor, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  estimated mean for measurement 1

- m2:

  estimated mean for measurement 2

- sd1:

  estimated standard deviation for measurement 1

- sd2:

  estimated standard deviation for measurement 2

- cor:

  estimated correlation between measurements

- n:

  sample size

## Value

Returns a 3-row matrix. The columns are:

- Estimate - estimated standardized mean difference

- adj Estimate - bias adjusted standardized mean difference estimate

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2008). “Confidence intervals for standardized linear
contrasts of means.” *Psychological Methods*, **13**(2), 99–109. ISSN
1939-1463,
[doi:10.1037/1082-989X.13.2.99](https://doi.org/10.1037/1082-989X.13.2.99)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.stdmean.ps(.05, 602.4, 705.6, 33.17, 51.08, .769, 8)
#>                             Estimate adj Estimate      SE      LL      UL
#> Unweighted standardizer:     -2.3963      -2.2185 0.65209 -3.6744 -1.1182
#> Measurement 1 standardizer:  -3.1112      -2.7656 0.91362 -4.9019 -1.3206
#> Measurement 2 standardizer:  -2.0204      -1.7959 0.59328 -3.1832 -0.8576

# Should return:
#                             Estimate adj Estimate      SE      LL      UL
# Unweighted standardizer:     -2.3963      -2.2185 0.65209 -3.6744 -1.1182
# Measurement 1 standardizer:  -3.1112      -2.7656 0.91362 -4.9019 -1.3206
# Measurement 2 standardizer:  -2.0204      -1.7959 0.59328 -3.1832 -0.8576

```
