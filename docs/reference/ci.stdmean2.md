# Confidence intervals for a 2-group standardized mean difference

Computes confidence intervals for a population standardized mean
difference. Unweighted, weighted, and single group variance
standardizers are used. The square root weighted variance standardizer
is recommended in 2-group nonexperimental designs with simple random
sampling. The square root unweighted variance standardizer is
recommended in 2-group experimental designs. The single group standard
deviation standardizer can be used with experimental or nonexperimental
designs. Equality of variances is not assumed.

For more details, see Section 2.4 of Bonett (2021, Volume 1)

## Usage

``` r
ci.stdmean2(alpha, m1, m2, sd1, sd2, n1, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  estimated mean for group 1

- m2:

  estimated mean for group 2

- sd1:

  estimated standard deviation for group 1

- sd2:

  estimated standard deviation for group 2

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a 4-row matrix. The columns are:

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
ci.stdmean2(.05, 20.9, 19.1, 3.85, 3.19, 50, 50)
#>                          Estimate adj Estimate      SE     LL     UL
#> Unweighted standardizer:   0.5091       0.5052 0.20539 0.1066 0.9117
#> Weighted standardizer:     0.5091       0.5052 0.20328 0.1107 0.9076
#> Group 1 standardizer:      0.4675       0.4603 0.19144 0.0923 0.8427
#> Group 2 standardizer:      0.5643       0.5556 0.23105 0.1114 1.0171

# Should return:
#                          Estimate adj Estimate      SE     LL     UL
# Unweighted standardizer:   0.5091       0.5052 0.20539 0.1066 0.9117
# Weighted standardizer:     0.5091       0.5052 0.20328 0.1107 0.9076
# Group 1 standardizer:      0.4675       0.4603 0.19144 0.0923 0.8427
# Group 2 standardizer:      0.5643       0.5556 0.23105 0.1114 1.0171

```
