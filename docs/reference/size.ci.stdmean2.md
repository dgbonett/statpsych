# Sample size for a 2-group standardized mean difference confidence interval

Computes the sample size per group required to estimate two types of
population standardized mean differences (unweighted standardizer and
single group standardizer) with desired confidence interval precision in
a 2-group design. Set the standardized mean difference planning value to
the largest value within a plausible range for a conservatively large
sample size. Set R = 1 for equal sample sizes. For unequal sample sizes,
this function assumes approximately equal population variances.

For more details, see Section 2.13 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.stdmean2(alpha, d, w, R)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- d:

  planning value of standardized mean difference

- w:

  desired confidence interval width

- R:

  n2/n1 ratio

## Value

Returns the required sample size per group for each standardizer

## References

Bonett DG (2009). “Estimating standardized linear contrasts of means
with desired precision.” *Psychological Methods*, **14**(1), 1–5. ISSN
1939-1463, [doi:10.1037/a0014270](https://doi.org/10.1037/a0014270) .
Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.stdmean2(.05, 1.0, .5, 1)
#>                             n1  n2
#> Unweighted standardizer:   139 139
#> Single group standardizer: 154 154

# Should return:
#                              n1  n2
# Unweighted standardizer:    139 139
# Single group standardizer:  154 154

size.ci.stdmean2(.05, 1.0, .5, 2)
#>                             n1  n2
#> Unweighted standardizer:   104 208
#> Single group standardizer: 116 232

# Should return:
#                              n1  n2
# Unweighted standardizer:    104 208
# Single group standardizer:  116 232
 
```
