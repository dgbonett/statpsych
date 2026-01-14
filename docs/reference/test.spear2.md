# Hypothesis test for a 2-group Spearman correlation difference

Computes a z test for a difference of population Spearman correlations
in a 2-group design. The test statistic uses a Bonett-Wright standard
error for each Spearman correlation. The hypothesis testing results
should be accompanied with a confidence interval for a difference in
population Spearman correlation values (see
[ci.spear2](https://dgbonett.github.io/statpsych/reference/ci.spear2.md)).

## Usage

``` r
test.spear2(cor1, cor2, n1, n2)
```

## Arguments

- cor1:

  estimated Spearman correlation for group 1

- cor2:

  estimated Spearman correlation for group 2

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of correlation difference

- z - z test statistic

- p - two-sided p-value

## References

Bonett DG, Wright TA (2000). “Sample size requirements for estimating
Pearson, Kendall and Spearman correlations.” *Psychometrika*, **65**(1),
23–28. ISSN 0033-3123,
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183) .

## Examples

``` r
test.spear2(.684, .437, 100, 125)
#>  Estimate      z       p
#>     0.247 2.4986 0.01247

# Should return:
# Estimate      z       p
#    0.247 2.4986 0.01247

```
