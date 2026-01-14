# Sample size for a paired-samples mean difference confidence interval

Computes the sample size required to estimate a difference in population
means with desired confidence interval precision in a paired-samples
design. Set the Pearson correlation planning value to the smallest value
within a plausible range for a conservatively large sample size. Set the
variance planning value to the largest value within a plausible range
for a conservatively large sample size.

For more details, see Section 4.26 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.mean.ps(alpha, var, cor, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  planning value of average variance of the two measurements

- cor:

  planning value of correlation between measurements

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.mean.ps(.05, 5.0, .5, 2.0)
#>  Sample size
#>           22

# Should return:
# Sample size
#          22
 
```
