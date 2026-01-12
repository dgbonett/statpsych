# Upper confidence limit of a variance

Computes an upper confidence limit for a population variance using an
estimated variance from a sample of size n in a prior study. The upper
limit can be used as a variance planning value in sample size functions
for desired power that require a planning value of the population
variance.

For more details, see Section 1.31 of Bonett (2021, Volume 1)

## Usage

``` r
ci.var.upper(alpha, var, n)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence (one-sided)

- var:

  estimated variance

- n:

  sample size

## Value

Returns an upper limit (UL) variance planning value

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.var.upper(.10, 1.45, 100)
#>        UL
#>  1.762447

# Should return:
#       UL
# 1.762447
 
```
