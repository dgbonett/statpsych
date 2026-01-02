# Sample size for a Spearman correlation confidence interval

Computes the sample size required to estimate a population Spearman
correlation with desired confidence interval precision. Set the
correlation planning value to the smallest absolute value within a
plausible range for a conservatively large sample size.

## Usage

``` r
size.ci.spear(alpha, cor, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor:

  planning value of Spearman correlation

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG, Wright TA (2000). “Sample size requirements for estimating
Pearson, Kendall and Spearman correlations.” *Psychometrika*, **65**(1),
23–28. ISSN 0033-3123,
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183) .

## Examples

``` r
size.ci.spear(.05, .362, .25)
#>  Sample size
#>          200

# Should return:
# Sample size
#         200
 
```
