# Sample size for a Pearson or partial correlation confidence interval

Computes the sample size required to estimate a population Pearson or
partial correlation with desired confidence interval precision. Set s =
0 for a Pearson correlation. Set the correlation planning value to the
smallest absolute value within a plausible range for a conservatively
large sample size.

For more details, see Section 1.24 of Bonett (2021, Volume 2)

## Usage

``` r
size.ci.cor(alpha, cor, s, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor:

  planning value of correlation

- s:

  number of control variables

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG, Wright TA (2000). “Sample size requirements for estimating
Pearson, Kendall and Spearman correlations.” *Psychometrika*, **65**(1),
23–28. ISSN 0033-3123,
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183) . Bonett DG
(2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.cor(.05, .3, 0, .2)
#>  Sample size
#>          320

# Should return:
# Sample size
#         320
 
```
