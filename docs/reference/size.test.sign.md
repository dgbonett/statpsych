# Sample size for a 1-group sign test

Computes the sample size required for a 1-group sign test with desired
power (see size.test.sign.ps for a paired-samples sign test). The Sign
test is a test of the null hypothesis that the population median is
equal to some specified value. This null hypothesis can also be
expressed in terms of the proportion of scores in the population that
are greater than the hypothesized population median value. Under the
null hypothesis, the population proportion is equal to .5. This function
requires a planning value of the population proportion.

For more details, see Section 1.29 of Bonett (2021, Volume 1)

## Usage

``` r
size.test.sign(alpha, pow, p)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- p:

  planning value of proportion

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.test.sign(.05, .90, .3)
#>  Sample size
#>           67

# Should return:
# Sample size
#          67
 
```
