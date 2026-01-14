# Sample size for a test of a single proportion

Computes the sample size required to test a population proportion with
desired power (using a correction for continuity) in a 1-group design.

For more details, see Section 1.15 of Bonett (2021, Volume 3)

## Usage

``` r
size.test.prop(alpha, pow, p, h)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- p:

  planning value of proportion

- h:

  null hypothesis value of proportion

## Value

Returns the required sample size

## References

Fleiss JL, Paik MC (2003). *Statistical Methods for Rates and
Proportions*, 3rd edition. Wiley.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.test.prop(.05, .9, .5, .3)
#>  Sample size
#>           65

# Should return:
# Sample size
#          65

```
