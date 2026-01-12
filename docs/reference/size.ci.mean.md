# Sample size for a mean confidence interval

Computes the sample size required to estimate a population mean with
desired confidence interval precision. Set the variance planning value
to the largest value within a plausible range for a conservatively large
sample size.

For more details, see Section 1.28 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.mean(alpha, var, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  planning value of response variable variance

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.mean(.05, 6.0, 1.5)
#>  Sample size
#>           43

# Should return:
# Sample size
#          43
 
```
