# Sample size for a paired-samples standardized mean difference confidence interval

Computes the sample size required to estimate two types of population
standardized mean differences (unweighted standardizer and single group
standardizer) with desired confidence interval precision in a
paired-samples design. Set the standardized mean difference planning
value to the largest value within a plausible range, and set the Pearson
correlation planning value to the smallest value within a plausible
range for a conservatively large sample size.

For more details, see Section 4.26 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.stdmean.ps(alpha, d, cor, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- d:

  planning value of standardized mean difference

- cor:

  planning value of correlation between measurements

- w:

  desired confidence interval width

## Value

Returns the required sample size for each standardizer

## References

Bonett DG (2009). “Estimating standardized linear contrasts of means
with desired precision.” *Psychological Methods*, **14**(1), 1–5. ISSN
1939-1463, [doi:10.1037/a0014270](https://doi.org/10.1037/a0014270) .
Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.stdmean.ps(.05, 1, .65, .6)
#>                            Sample size
#> Unweighted standardizer:            46
#> Single group standardizer:          52

# Should return:
#                            Sample Size
# Unweighted standardizer:            46
# Single group standardizer:          52
 
```
