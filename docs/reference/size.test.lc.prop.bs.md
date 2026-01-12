# Sample size for a test of between-subjects proportion linear contrast

Computes the sample size in each group (assuming equal sample sizes)
required to test a linear contrast of population proportions with
desired power in a between-subjects design. The planning value for the
effect size (linear contrast of proportions) could be set equal to the
linear contrast of proportion planning values or it could be set equal
to a minimally interesting effect size. For a conservatively large
sample size, set the proportion planning values to .5 and set the effect
size to a minimally interesting value.

For more details, see Section 2.23 of Bonett (2021, Volume 3)

## Usage

``` r
size.test.lc.prop.bs(alpha, pow, p, es, v)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- p:

  vector of proportion planning values

- es:

  planning value of proportion linear contrast

- v:

  vector of between-subjects contrast coefficients

## Value

Returns the required sample size for each group

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
p <- c(.35, .35, .20)
v <- c(.5, .5, -1)
size.test.lc.prop.bs(.05, .9, p, .15, v)
#>  Sample size per group
#>                    128

# Should return:
# Sample size per group
#                   128

```
