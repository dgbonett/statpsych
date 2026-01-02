# Confidence interval for a 2-group Pearson correlation difference

Computes a confidence interval for a difference in population Pearson
correlations in a 2-group design. A bias adjustment is used to reduce
the bias of each Fisher transformed correlation.

For more details, see Section 2.17 of Bonett (2021, Volume 2)

## Usage

``` r
ci.cor2(alpha, cor1, cor2, n1, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor1:

  estimated Pearson correlation for group 1

- cor2:

  estimated Pearson correlation for group 2

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

Zou GY (2007). “Toward using confidence intervals to compare
correlations.” *Psychological Methods*, **12**(4), 399–413. ISSN
1939-1463,
[doi:10.1037/1082-989X.12.4.399](https://doi.org/10.1037/1082-989X.12.4.399)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.cor2(.05, .64, .31, 200, 200)
#>  Estimate      SE     LL     UL
#>      0.33 0.07692 0.1797 0.4814

# Should return:
# Estimate      SE      LL     UL
#     0.33 0.07692  0.1797 0.4814
 
```
