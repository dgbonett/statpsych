# Confidence interval for a 2-group Spearman correlation difference

Computes a confidence interval for a difference of population Spearman
correlations in a 2-group design. This function is not appropriate for
ordered categorical variables.

## Usage

``` r
ci.spear2(alpha, cor1, cor2, n1, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

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

- Estimate - estimated correlation difference

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Wright TA (2000). “Sample size requirements for estimating
Pearson, Kendall and Spearman correlations.” *Psychometrika*, **65**(1),
23–28. ISSN 0033-3123,
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183) .

Zou GY (2007). “Toward using confidence intervals to compare
correlations.” *Psychological Methods*, **12**(4), 399–413. ISSN
1939-1463,
[doi:10.1037/1082-989X.12.4.399](https://doi.org/10.1037/1082-989X.12.4.399)
.

## Examples

``` r
ci.spear2(.05, .54, .48, 180, 200)
#>  Estimate      SE      LL     UL
#>      0.06 0.08125 -0.1004 0.2185

# Should return:
# Estimate      SE      LL     UL
#     0.06 0.08125 -0.1004 0.2185     
 
```
