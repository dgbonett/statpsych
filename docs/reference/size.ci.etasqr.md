# Sample size for an eta-squared confidence interval

Computes the sample size required to estimate an eta-squared coefficient
in a one-way ANOVA with desired confidence interval precision. Set the
planning value of eta-squared to about 1/3 for a conservatively large
sample size.

For more details, see Section 3.24 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.etasqr(alpha, etasqr, groups, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- etasqr:

  planning value of eta-squared

- groups:

  number of groups

- w:

  desired confidence interval width

## Value

Returns the required sample size for each group

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.etasqr(.05, .2, 3, .15)
#>  Sample size per group
#>                    103

# Should return:
# Sample size per group 
#                   103
 
```
