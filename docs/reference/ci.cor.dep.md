# Confidence interval for a difference in dependent Pearson correlations

Computes a confidence interval for a difference in population Pearson
correlations that are estimated from the same sample and have one
variable in common. A bias adjustment is used to reduce the bias of each
Fisher transformed correlation. An approximate standard error is
recovered from the confidence interval.

For more details, see Section 2.17 of Bonett (2021, Volume 2)

## Usage

``` r
ci.cor.dep(alpha, cor1, cor2, cor12, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor1:

  estimated Pearson correlation between y and x1

- cor2:

  estimated Pearson correlation between y and x2

- cor12:

  estimated Pearson correlation between x1 and x2

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated correlation difference

- SE - recovered standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Zou GY (2007). “Toward using confidence intervals to compare
correlations.” *Psychological Methods*, **12**(4), 399–413. ISSN
1939-1463,
[doi:10.1037/1082-989X.12.4.399](https://doi.org/10.1037/1082-989X.12.4.399)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.cor.dep(.05, .396, .179, .088, 166)
#>  Estimate     SE     LL     UL
#>     0.217 0.1027 0.0132 0.4158

# Should return:
# Estimate     SE     LL     UL
#    0.217 0.1027 0.0132 0.4158
 
```
