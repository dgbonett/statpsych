# Sample size for a paired-samples proportion ratio confidence interval

Computes the sample size required to estimate a ratio of population
proportions with desired confidence interval precision in a
paired-samples design. Set the phi correlation planning value to the
smallest value within a plausible range for a conservatively large
sample size.

For more details, see Section 3.12 of Bonett (2021, Volume 3)

## Usage

``` r
size.ci.ratio.prop.ps(alpha, p1, p2, phi, r)
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

- r:

  desired upper to lower confidence interval endpoint ratio

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.ratio.prop.ps(.05, .4, .2, .7, 2)
#>  Sample size
#>           67

# Should return:
# Sample size
#          67

```
