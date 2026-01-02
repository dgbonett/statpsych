# Confidence interval for a 2-group correlation difference

Computes a 100(1 - alpha)% confidence interval for a difference in
population correlations in a 2-group design. The correlations can be
Pearson, Spearman, partial, semipartial, or point-biserial correlations.
The correlations could also be correlations between two latent factors.
The function requires a point estimate and a 100(1 - alpha)% confidence
interval for each correlation as input. The confidence intervals for
each correlation can be obtained using the ci.fisher function.

For more details, see Section 2.17 of Bonett (2021, Volume 2)

## Usage

``` r
ci.cor2.gen(cor1, ll1, ul1, cor2, ll2, ul2)
```

## Arguments

- cor1:

  estimated correlation for group 1

- ll1:

  lower limit for group 1 correlation

- ul1:

  upper limit for group 1 correlation

- cor2:

  estimated correlation for group 2

- ll2:

  lower limit for group 2 correlation

- ul2:

  upper limit for group 2 correlation

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated correlation difference

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
ci.cor2.gen(.64, .55, .71, .31, .18, .43)
#>  Estimate   LL     UL
#>      0.33 0.18 0.4776

# Should return:
# Estimate    LL     UL
#      0.33 0.18 0.4776
 
```
