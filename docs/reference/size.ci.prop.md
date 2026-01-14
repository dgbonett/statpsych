# Sample size for a proportion confidence interval

Computes the sample size required to estimate a population proportion
with desired confidence interval precision. Set the proportion planning
value to .5 for a conservatively large sample size.

For more details, see Section 1.14 of Bonett (2021, Volume 3)

## Usage

``` r
size.ci.prop(alpha, p, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- p:

  planning value of proportion

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.prop(.05, .5, .1)
#>  Sample size
#>          385

# Should return:
# Sample size
#         385

```
