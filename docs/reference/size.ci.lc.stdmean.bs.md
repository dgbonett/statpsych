# Sample size for a between-subjects standardized linear contrast of means confidence interval

Computes the sample size per group (assuming equal sample sizes)
required to estimate two types of standardized linear contrasts of
population means (unweighted average standardizer and single group
standardizer) with desired confidence interval precision in a
between-subjects design. Set the standardized linear contrast of means
to the largest value within a plausible range for a conservatively large
sample size.

For more details, see Section 3.24 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.lc.stdmean.bs(alpha, d, w, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- d:

  planning value of standardized linear contrast of means

- w:

  desired confidence interval width

- v:

  vector of between-subjects contrast coefficients

## Value

Returns the required sample size per group for each standardizer

## References

Bonett DG (2009). “Estimating standardized linear contrasts of means
with desired precision.” *Psychological Methods*, **14**(1), 1–5. ISSN
1939-1463, [doi:10.1037/a0014270](https://doi.org/10.1037/a0014270) .
Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
v <- c(.5, .5, -1)
size.ci.lc.stdmean.bs(.05, .8, .6, v)
#>                            Sample size per group
#> Unweighted standardizer:                      69
#> Single group standardizer:                    78

# Should return:
#                            Sample size per group
# Unweighted standardizer:                      69
# Single group standardizer:                    78
 
```
