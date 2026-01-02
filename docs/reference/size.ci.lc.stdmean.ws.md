# Sample size for a within-subjects standardized linear contrast of means confidence interval

Computes the sample size required to estimate two types of standardized
linear contrasts of population means (unweighted standardizer and single
level standardizer) with desired confidence interval precision in a
within-subjects design. For a conservatively large sample size, set the
standardized linear contrast of means planning value to the largest
value within a plausible range, and set the Pearson correlation planning
value to the smallest value within a plausible range.

For more details, see Section 4.26 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.lc.stdmean.ws(alpha, d, cor, w, q)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- d:

  planning value of standardized linear contrast

- cor:

  planning value of average correlation between measurements

- w:

  desired confidence interval width

- q:

  vector of within-subjects contrast coefficients

## Value

Returns the required sample size for each standardizer

## References

Bonett DG (2009). “Estimating standardized linear contrasts of means
with desired precision.” *Psychological Methods*, **14**(1), 1–5. ISSN
1939-1463, [doi:10.1037/a0014270](https://doi.org/10.1037/a0014270) .
Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
q <- c(1/3, 1/3, 1/3, -1)
size.ci.lc.stdmean.ws(.05, .5, .7, .4, q)
#>                            Sample size
#> Unweighted standardizer:            46
#> Single level standardizer:          51

# Should return:
#                            Sample size
# Unweighted standardizer:            46
# Single level standardizer:          51
 
```
