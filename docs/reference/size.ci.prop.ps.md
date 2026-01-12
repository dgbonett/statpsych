# Sample size for a paired-sample proportion difference confidence interval

Computes the sample size required to estimate a population proportion
difference with desired confidence interval precision in a
paired-samples design. Set the proportion planning values to .5 for a
conservatively large sample size. Set the phi correlation planning value
to the smallest value within a plausible range for a conservatively
large sample size.

For more details, see Section 3.12 of Bonett (2021, Volume 3)

## Usage

``` r
size.ci.prop.ps(alpha, p1, p2, phi, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- p1:

  planning value of proportion for measurement 1

- p2:

  planning value of proportion for measurement 2

- phi:

  planning value of phi correlation

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.prop.ps(.05, .25, .35, .6, .1)
#>  Sample size
#>          257

# Should return:
# Sample size
#         257

```
