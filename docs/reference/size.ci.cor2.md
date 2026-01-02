# Sample size for a 2-group Pearson correlation difference confidence interval

Computes the sample size required to estimate a difference in population
Pearson or partial correlations with desired confidence interval
precision in a 2-group design. Set the correlation planning values to
the smallest absolute values within their plausible ranges for a
conservatively large sample size.

## Usage

``` r
size.ci.cor2(alpha, cor1, cor2, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor1:

  correlation planning value for group 1

- cor2:

  correlation planning value for group 2

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
size.ci.cor2(.05, .8, .5, .2)
#>  Sample size per group
#>                    271

# Should return:
# Sample size per group
#                   271
 
```
