# Sample size for a test of a Pearson or partial correlation

Computes the sample size required to test a population Pearson or a
partial correlation with desired power. Set s = 0 for a Pearson
correlation.

For more details, see Section 1.25 of Bonett (2021, Volume 2)

## Usage

``` r
size.test.cor(alpha, pow, cor, s, h)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- cor:

  planning value of correlation

- s:

  number of control variables

- h:

  null hypothesis value of correlation

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.test.cor(.05, .95, -.50, 0, 0)
#>  Sample size
#>           47

# Should return:
# Sample size
#          47

size.test.cor(.05, .90, .4, 2, 0)
#>  Sample size
#>           64

# Should return:
# Sample size
#          64
 
```
